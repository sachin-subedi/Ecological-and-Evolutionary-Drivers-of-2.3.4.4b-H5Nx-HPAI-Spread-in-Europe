import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import pickle
from pathlib import Path

ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR"
]
CLC_RASTER    = "U2018_CLC2018_V2020_20u1.tif"
SHAPEFILE = "ne_110m_admin_0_countries.shp"
LOOKUP_CSV    = "cluster_lookup_HC.csv"
CELL_AREA_KM2 = 0.01
PICKLE_PATH   = "filtered_stats.pkl"
MIN_AREA      = 0.00001   # substitute for true zeros

AI_CLC_CODES = [
    10,12,13,14,18,21,22,23,24,25,27,28,29,30,31,32,35,36,37,38,39,40,41,42,43,44
]

CAT_MAP = {
    10:"AgriculturalLand", 12:"AgriculturalLand", 13:"AgriculturalLand",
    14:"AgriculturalLand", 18:"AgriculturalLand", 21:"AgriculturalLand",
    22:"AgriculturalLand",
    23:"ForestAreas", 24:"ForestAreas", 25:"ForestAreas", 27:"ForestAreas",
    28:"ForestAreas", 30:"ForestAreas", 31:"ForestAreas", 32:"ForestAreas",
    37:"Wetlands", 38:"Wetlands", 39:"Wetlands",
    40:"WaterBodies", 41:"WaterBodies", 42:"WaterBodies",
    43:"WaterBodies", 44:"WaterBodies"
}

ALL_CATEGORIES = ["AgriculturalLand", "ForestAreas", "Wetlands", "WaterBodies"]
ALL_CLUSTERS   = [
    "HC_Cluster_1_Alpine", "HC_Cluster_1_Atlantic", "HC_Cluster_1_Continental",
    "HC_Cluster_2_Alpine", "HC_Cluster_2_Continental", "HC_Cluster_2_Mediterranean",
    "HC_Cluster_2_Pannonian", "HC_Cluster_3_Alpine", "HC_Cluster_3_Boreal",
    "HC_Cluster_4_Mediterranean"
]

gdf_raw = gpd.read_file(SHAPEFILE)

gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Norway", "ISO_A3"] = "NOR"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "France",  "ISO_A3"] = "FRA"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Kosovo",  "ISO_A3"] = "XKX"

countries = (
    gdf_raw
    .loc[lambda d: d["ISO_A3"].isin(ISO38), ["ISO_A3", "SOVEREIGNT", "geometry"]]
    .rename(columns={"ISO_A3": "ISO3", "SOVEREIGNT": "Country"})
    .reset_index(drop=True)
)

if countries.crs is None:
    countries = countries.set_crs(epsg=4326)
else:
    countries = countries.to_crs(epsg=4326)
countries = countries.to_crs(epsg=3035)
countries = countries[countries["Country"] != "Montenegro"]

print(f"Countries loaded: {len(countries)}")
print(countries[["ISO3", "Country"]])

cluster_lookup = (
    pd.read_csv(LOOKUP_CSV)
    .rename(columns=str.strip)
    .rename(columns={"GeoCluster": "cluster"})
)
countries = countries.merge(cluster_lookup, on="Country", how="left")

if countries["cluster"].isna().any():
    missing = countries.loc[countries["cluster"].isna(), "Country"].tolist()
    raise ValueError(f"No GeoCluster found for: {', '.join(missing)}")

if Path(PICKLE_PATH).exists():
    print("Loading cached zonal stats...")
    with open(PICKLE_PATH, "rb") as f:
        filtered_stats = pickle.load(f)
else:
    print("Running zonal stats (this may take a while)...")
    raw_stats = zonal_stats(
        countries["geometry"],
        CLC_RASTER,
        categorical=True,
        nodata=0
    )
    filtered_stats = [
        {code: cnt for code, cnt in d.items() if code in AI_CLC_CODES}
        for d in raw_stats
    ]
    with open(PICKLE_PATH, "wb") as f:
        pickle.dump(filtered_stats, f)

country_long = pd.DataFrame(
    [
        {"ISO3": iso, "cluster": clu, "CLC_code": code, "count": cnt}
        for iso, clu, d in zip(countries.ISO3, countries.cluster, filtered_stats)
        for code, cnt in d.items()
    ]
)
country_long["category"] = country_long["CLC_code"].map(CAT_MAP)

cluster_cat_cnt = (
    country_long
    .groupby(["cluster", "category"], as_index=False)["count"]
    .sum()
)
cluster_cat_cnt["area_km2"] = cluster_cat_cnt["count"] * CELL_AREA_KM2

full_index = pd.MultiIndex.from_product(
    [ALL_CLUSTERS, ALL_CATEGORIES], names=["cluster", "category"]
)
cluster_cat_cnt = (
    cluster_cat_cnt
    .set_index(["cluster", "category"])
    .reindex(full_index)
    .reset_index()
)

cluster_cat_cnt["area_km2"] = cluster_cat_cnt["area_km2"].fillna(MIN_AREA)
cluster_cat_cnt = cluster_cat_cnt.drop(columns=["count"], errors="ignore")

wide_area = (
    cluster_cat_cnt
    .pivot(index="cluster", columns="category", values="area_km2")
    .fillna(MIN_AREA)
    .sort_index()
)

print("\nArea by land-cover category (km²):")
print(wide_area.round(5))

for cat in ALL_CATEGORIES:
    area_df = (
        cluster_cat_cnt.loc[cluster_cat_cnt["category"] == cat, ["cluster", "area_km2"]]
        .rename(columns={"cluster": "GeoCluster", "area_km2": f"{cat}_areas"})
        .sort_values("GeoCluster")
    )
    area_df.to_csv(f"{cat}_areas.tsv", sep="\t", index=False)

for cat in ALL_CATEGORIES:
    print(f" • {cat}_areas.tsv")
