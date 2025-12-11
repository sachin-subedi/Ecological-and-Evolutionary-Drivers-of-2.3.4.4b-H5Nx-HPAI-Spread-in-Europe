import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import xarray as xr, rioxarray

GRIB_FILE    = "data.grib"
SHAPEFILE    = "/home/ss11645/Landcover/shapes/ne_110m_admin_0_countries.shp"
LAI_LV_TIF   = "lai_lv_mean_2016_2025.tif"
LAI_HV_TIF   = "lai_hv_mean_2016_2025.tif"

ds = xr.open_dataset(GRIB_FILE, engine="cfgrib")

lai_lv_mean = ds['lai_lv'].mean("time", skipna=True).rio.write_crs("EPSG:4326")
lai_hv_mean = ds['lai_hv'].mean("time", skipna=True).rio.write_crs("EPSG:4326")

lai_lv_mean.rio.to_raster(LAI_LV_TIF, compress="LZW")
lai_hv_mean.rio.to_raster(LAI_HV_TIF, compress="LZW")
print("✔️  GeoTIFFs written")

gdf = gpd.read_file(SHAPEFILE)

ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR","BLR"
]
gdf = gdf.set_crs("EPSG:4326")

gdf = gdf[gdf["ISO_A3"].isin(ISO38)].copy()
print(f"Countries retained: {len(gdf)}")

cluster_map = {
    "One":   ["ALB","BGR","CYP","GRC","XKX","MKD","SRB","ROU"],
    "Two":   ["HUN","SVK","POL","UKR","MDA","BLR"],
    "Three": ["DNK","FRA","BEL","DEU","ISL","IRL","LUX","NLD","PRT","ESP","CHE","GBR"],
    "Four":  ["EST","FIN","LVA","LTU","NOR","SWE"],
    "Five":  ["AUT","MNE","BIH","HRV","CZE","ITA","SVN"],
}

iso_to_cluster = {iso: z for z, lst in cluster_map.items() for iso in lst}
gdf["GeoCluster"] = gdf["ISO_A3"].map(iso_to_cluster)
print(gdf[["ISO_A3","GeoCluster"]].head())

# LAI_LV
stats_lv = zonal_stats(
    gdf, LAI_LV_TIF, stats=["mean"], geojson_out=True, nodata=-999
)
lv_df = gpd.GeoDataFrame.from_features(stats_lv)[["ISO_A3","mean"]].rename(columns={"mean":"LAI_LV"})

# LAI_HV
stats_hv = zonal_stats(
    gdf, LAI_HV_TIF, stats=["mean"], geojson_out=True, nodata=-999
)
hv_df = gpd.GeoDataFrame.from_features(stats_hv)[["ISO_A3","mean"]].rename(columns={"mean":"LAI_HV"})

# merge back
df = gdf[["ISO_A3","GeoCluster"]].merge(lv_df,on="ISO_A3").merge(hv_df,on="ISO_A3")
print(df.head())


cluster_means = (
    df.groupby("GeoCluster")[["LAI_LV","LAI_HV"]]
      .mean()
      .round(3)
      .sort_index()
)
print("\n=== Mean LAI 2016–2025 by GeoCluster ===")
print(cluster_means)

cluster_means.to_csv("GeoCluster_LAI_mean_2016_2025.csv")

cluster_means = (
    df.groupby("GeoCluster")[["LAI_LV", "LAI_HV"]]
      .mean()
      .round(3)
      .sort_index()
)

cluster_means[["LAI_LV"]].to_csv("GeoCluster_LAI_LV.tsv", sep="\t")

cluster_means[["LAI_HV"]].to_csv("GeoCluster_LAI_HV.tsv", sep="\t")

