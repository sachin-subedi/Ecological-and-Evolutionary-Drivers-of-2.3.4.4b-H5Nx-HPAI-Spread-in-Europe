import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import xarray as xr, rioxarray

GRIB_FILE  = "data.grib"
SHAPEFILE  = "ne_110m_admin_0_countries.shp"
LAI_LV_TIF = "lai_lv_mean_2016_2025.tif"
LAI_HV_TIF = "lai_hv_mean_2016_2025.tif"

ds = xr.open_dataset(GRIB_FILE, engine="cfgrib")
lai_lv_mean = ds['lai_lv'].mean("time", skipna=True).rio.write_crs("EPSG:4326")
lai_hv_mean = ds['lai_hv'].mean("time", skipna=True).rio.write_crs("EPSG:4326")
lai_lv_mean.rio.to_raster(LAI_LV_TIF, compress="LZW")
lai_hv_mean.rio.to_raster(LAI_HV_TIF, compress="LZW")
print("✔️  GeoTIFFs written")

ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR","BLR"
]

hc_cluster_map = {
    "HC_Cluster_1_Alpine":          ["CHE"],
    "HC_Cluster_1_Atlantic":        ["BEL", "FRA", "NLD", "GBR", "ISL", "IRL"],
    "HC_Cluster_1_Continental":     ["DEU", "LUX"],
    "HC_Cluster_2_Alpine":          ["AUT", "SVK", "BIH"],
    "HC_Cluster_2_Continental":     ["BGR", "HRV", "CZE", "XKX", "MDA", "POL",
                                     "ROU", "SVN", "DNK", "MKD", "SRB", "UKR"],
    "HC_Cluster_2_Mediterranean":   ["ALB", "GRC", "ITA", "CYP"],
    "HC_Cluster_2_Pannonian":       ["HUN"],
    "HC_Cluster_3_Alpine":          ["NOR"],
    "HC_Cluster_3_Boreal":          ["EST", "FIN", "LVA", "LTU", "SWE"],
    "HC_Cluster_4_Mediterranean":   ["ESP", "PRT"],
}

iso_to_cluster = {iso: cluster for cluster, members in hc_cluster_map.items() for iso in members}

gdf_raw = gpd.read_file(SHAPEFILE)

gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Norway", "ISO_A3"] = "NOR"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "France", "ISO_A3"] = "FRA"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Kosovo", "ISO_A3"] = "XKX"

gdf = (
    gdf_raw
    .loc[lambda d: d["ISO_A3"].isin(ISO38)]
    .copy()
    .set_crs("EPSG:4326")
)

gdf["HC_Cluster"] = gdf["ISO_A3"].map(iso_to_cluster)

unmatched = gdf[gdf["HC_Cluster"].isna()]["ISO_A3"].tolist()
if unmatched:
    print(f"⚠️  No HC cluster assigned for: {unmatched}")

print(f"Countries retained: {len(gdf)}")
print(gdf[["ISO_A3", "HC_Cluster"]])

stats_lv = zonal_stats(gdf, LAI_LV_TIF, stats=["mean"], geojson_out=True, nodata=-999)
lv_df = gpd.GeoDataFrame.from_features(stats_lv)[["ISO_A3", "mean"]].rename(columns={"mean": "LAI_LV"})

stats_hv = zonal_stats(gdf, LAI_HV_TIF, stats=["mean"], geojson_out=True, nodata=-999)
hv_df = gpd.GeoDataFrame.from_features(stats_hv)[["ISO_A3", "mean"]].rename(columns={"mean": "LAI_HV"})

df = (
    gdf[["ISO_A3", "HC_Cluster"]]
    .merge(lv_df, on="ISO_A3")
    .merge(hv_df, on="ISO_A3")
)
print(df.head())

cluster_means = (
    df.groupby("HC_Cluster")[["LAI_LV", "LAI_HV"]]
      .mean()
      .round(3)
      .sort_index()
)

print("\n=== Mean LAI 2016–2025 by HC_Cluster ===")
print(cluster_means)

cluster_means.to_csv("HC_GeoCluster_LAI_mean_2016_2025.csv")
cluster_means[["LAI_LV"]].to_csv("HC_GeoCluster_LAI_LV.tsv", sep="\t")
cluster_means[["LAI_HV"]].to_csv("HC_GeoCluster_LAI_HV.tsv", sep="\t")
