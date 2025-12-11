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
CLC_RASTER = "/home/ss11645/Landcover/U2018_CLC2018_V2020_20u1.tif"
SHAPEFILE = "/home/ss11645/Landcover/shapes/ne_110m_admin_0_countries.shp"
LOOKUP_CSV = "cluster_lookup_auto.csv"
CELL_AREA_KM2 = 0.01  # 100 m × 100 m
PICKLE_PATH = "filtered_stats.pkl"

AI_CLC_CODES = [
    10,12,13,14,18,21,22,23,24,25,27,28,29,30,31,32,35,36,37,38,39,40,41,42,43,44
]

CAT_MAP = {
    # Agricultural land
    10:"AgriculturalLand",
    12:"AgriculturalLand",
    13:"AgriculturalLand",
    14:"AgriculturalLand",
    18:"AgriculturalLand",
    21:"AgriculturalLand",
    22:"AgriculturalLand",
    # Forest
    23:"ForestAreas",
    24:"ForestAreas",
    25:"ForestAreas",
    27:"ForestAreas",
    28:"ForestAreas",
    30:"ForestAreas",
    31:"ForestAreas",
    32:"ForestAreas",
    # Wetlands
    37:"Wetlands",
    38:"Wetlands",
    39:"Wetlands",
    # Water
    40:"WaterBodies",
    41:"WaterBodies",
    42:"WaterBodies",
    43:"WaterBodies",
    44:"WaterBodies"
}


countries = (
    gpd.read_file(SHAPEFILE)
      .set_crs(epsg=4326).to_crs(epsg=3035)  # match CLC raster CRS
      .loc[lambda d: d["ISO_A3_EH"].isin(ISO38), ["ISO_A3_EH", "ADMIN", "geometry"]]
      .rename(columns={"ISO_A3_EH": "ISO3", "ADMIN": "Country"})
      .reset_index(drop=True)
)
countries = countries[countries["Country"] != "Montenegro"]

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
    with open(PICKLE_PATH, "rb") as f:
        filtered_stats = pickle.load(f)
else:
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

totals = (
    cluster_cat_cnt.groupby("cluster", as_index=False)["count"]
      .sum()
      .rename(columns={"count": "total"})
)
cluster_cat_prop = (
    cluster_cat_cnt
      .merge(totals, on="cluster")
      .assign(prop=lambda d: d["count"] / d["total"])
      .drop(columns=["count", "total"])
)


cluster_cat_cnt[["cluster", "category", "area_km2"]].to_csv(
    "cluster_category_area.csv", index=False
)
cluster_cat_prop[["cluster", "category", "prop"]].to_csv(
    "cluster_category_proportion.csv", index=False
)

print("✅ Saved:")
print(" • cluster_category_area.csv")
print(" • cluster_category_proportion.csv")

wide_area = (
    cluster_cat_cnt
      .pivot(index="cluster", columns="category", values="area_km2")
      .fillna(0)
      .sort_index()
)
wide_prop = (
    cluster_cat_prop
      .pivot(index="cluster", columns="category", values="prop")
      .fillna(0)
      .sort_index()
)

print("\nArea by land-cover category (km²):")
print(wide_area.round(1))

print("\nProportion of each GeoCluster occupied by category:")
print((wide_prop * 100).round(2).astype(str) + " %")


for cat in ["AgriculturalLand", "ForestAreas", "Wetlands", "WaterBodies"]:
    area_df = (
        cluster_cat_cnt.loc[cluster_cat_cnt["category"] == cat, ["cluster", "area_km2"]]
          .rename(columns={
              "cluster": "GeoCluster",
              "area_km2": f"{cat}_areas"
          })
          .sort_values("GeoCluster")
    )
    area_df.to_csv(f"{cat}_areas.tsv", sep="\t", index=False)

    prop_df = (
        cluster_cat_prop.loc[cluster_cat_prop["category"] == cat, ["cluster", "prop"]]
          .rename(columns={
              "cluster": "GeoCluster",
              "prop": f"{cat}_prop"
          })
          .sort_values("GeoCluster")
    )
    prop_df.to_csv(f"{cat}_prop.tsv", sep="\t", index=False)

print("✅ Extra TSV files saved:")
print(" • AgriculturalLand_areas.tsv • AgriculturalLand_prop.tsv")
print(" • ForestAreas_areas.tsv • ForestAreas_prop.tsv")
print(" • Wetlands_areas.tsv • Wetlands_prop.tsv")
print(" • WaterBodies_areas.tsv • WaterBodies_prop.tsv")
