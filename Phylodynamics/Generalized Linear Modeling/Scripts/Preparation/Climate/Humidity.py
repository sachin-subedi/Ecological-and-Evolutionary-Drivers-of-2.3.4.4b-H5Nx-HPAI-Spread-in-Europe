import xarray as xr
import rioxarray
import geopandas as gpd
import regionmask
import pandas as pd
import numpy as np

ds = xr.open_dataset(
    "hu_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc",
    engine="netcdf4",
    chunks={"time": 365}
)
ds = ds.rio.write_crs(4326)
hu = ds["hu"]
print(ds)

ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR","BLR"
]

cluster_map = {
    "HC_Cluster_1_Alpine":        ["CHE"],
    "HC_Cluster_1_Atlantic":      ["BEL", "FRA", "NLD", "GBR", "ISL", "BLR"],
    "HC_Cluster_1_Continental":   ["DNK", "DEU", "LUX"],
    "HC_Cluster_2_Alpine":        ["AUT", "SVK", "BIH"],
    "HC_Cluster_2_Continental":   ["BGR", "HRV", "CZE", "XKX", "MDA", "POL", "ROU", "SVN"],
    "HC_Cluster_2_Mediterranean": ["ALB", "GRC", "ITA"],
    "HC_Cluster_2_Pannonian":     ["HUN"],
    "HC_Cluster_3_Alpine":        ["NOR"],
    "HC_Cluster_3_Boreal":        ["EST", "FIN", "LVA", "LTU", "SWE"],
    "HC_Cluster_4_Mediterranean": ["ESP", "PRT"]
}
iso_to_cluster = {iso: cluster for cluster, members in cluster_map.items() for iso in members}

SHAPEFILE = "ne_110m_admin_0_countries.shp"

gdf_raw = gpd.read_file(SHAPEFILE)

gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Norway", "ISO_A3"] = "NOR"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "France", "ISO_A3"] = "FRA"
gdf_raw.loc[gdf_raw["SOVEREIGNT"] == "Kosovo", "ISO_A3"] = "XKX"

gdf_country = (
    gdf_raw
    .loc[lambda d: d["ISO_A3"].isin(ISO38)]
    .assign(GeoCluster=lambda d: d["ISO_A3"].map(iso_to_cluster))
    .dropna(subset=["GeoCluster"])
)

print(gdf_country[["SOVEREIGNT", "ISO_A3", "GeoCluster"]])

if gdf_country.crs is None:
    gdf_country = gdf_country.set_crs("EPSG:4326")
else:
    gdf_country = gdf_country.to_crs("EPSG:4326")

country_cluster_lookup = gdf_country[["SOVEREIGNT", "ISO_A3", "GeoCluster"]].drop_duplicates()

gdf = gdf_country.dissolve(by="GeoCluster").reset_index()
gdf = gdf.reset_index(drop=True)
print(gdf[["GeoCluster"]])

mask2d = regionmask.mask_geopandas(gdf, hu.longitude, hu.latitude)
print("Unique mask values:", np.unique(mask2d.values[~np.isnan(mask2d.values)]))

cluster_means_list = []
for row_idx, row in gdf.iterrows():
    mean_hu = hu.where(mask2d == row_idx).mean(dim=("latitude", "longitude"))
    cluster_means_list.append(mean_hu.assign_coords(GeoCluster=row.GeoCluster))

cluster_means = xr.concat(cluster_means_list, dim="GeoCluster")

out_df = cluster_means.to_dataframe(name="hu_pct").reset_index()

out_with_country = country_cluster_lookup.merge(out_df, on="GeoCluster", how="left")

out_with_country.to_csv("GeoCluster_humidity_with_countries.tsv", sep="\t", index=False)

df_daily = (
    out_with_country
    .rename(columns={"SOVEREIGNT": "Country", "hu_pct": "humidity"})
    .loc[:, ["Country", "GeoCluster", "time", "humidity"]]
)

df_daily["time"] = pd.to_datetime(df_daily["time"])
df_daily = df_daily[
    (df_daily["time"] >= "2016-01-01") &
    (df_daily["time"] <= "2025-05-01")
]
df_daily["year"] = df_daily["time"].dt.year
df_daily["month"] = df_daily["time"].dt.month

monthly_country = (
    df_daily
    .groupby(["Country", "GeoCluster", "year", "month"], as_index=False)
    .agg(humidity=("humidity", "mean"))
)

avg_country = (
    monthly_country
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(humidity=("humidity", "mean"))
)

season_country = (
    monthly_country
    .groupby(["Country", "GeoCluster", "year"])
    .agg(humidity_sd=("humidity", "std"))
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(humidity_seasonality=("humidity_sd", "mean"))
)

avg_geo = (
    avg_country
    .groupby("GeoCluster", as_index=False)
    .agg(humidity=("humidity", "mean"))
)

print("\nCountry-level humidity mean:")
print(avg_country)

print("\nCountry-level humidity seasonality:")
print(season_country)

print("\nGeoCluster-level humidity mean:")
print(avg_geo)

avg_geo[["GeoCluster", "humidity"]].to_csv(
    "GeoCluster_humidity.tsv",
    sep="\t",
    index=False
)
print(avg_geo)
