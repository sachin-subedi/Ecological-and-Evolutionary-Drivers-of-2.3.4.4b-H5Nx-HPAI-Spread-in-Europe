import xarray as xr
import rioxarray
import geopandas as gpd
import regionmask
import pandas as pd
import xarray as xr

ds = xr.open_dataset(
    "tg_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc",
    engine="netcdf4",
    chunks={"time": 365}
)

ds = ds.rio.write_crs(4326)
tg = ds["tg"]
print(ds)

ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR","BLR"
]

cluster_map = {
    "HC_Cluster_1_Alpine":   ["CHE"],
    "HC_Cluster_1_Atlantic":   ["BEL", "FRA", "NLD", "GBR", "ISL", "BLR"],
    "HC_Cluster_1_Continental": ["DNK", "DEU", "LUX"],
    "HC_Cluster_2_Alpine":  ["AUT", "SVK", "BIH"],
    "HC_Cluster_2_Continental":  ["BGR", "HRV", "CZE", "XKX", "MDA", "POL", "ROU", "SVN"],
    "HC_Cluster_2_Mediterranean":   ["ALB", "GRC", "ITA"],
    "HC_Cluster_2_Pannonian":   ["HUN"],
    "HC_Cluster_3_Alpine": ["NOR"],
    "HC_Cluster_3_Boreal":  ["EST", "FIN", "LVA", "LTU", "SWE"],
    "HC_Cluster_4_Mediterranean":  ["ESP", "PRT"]
}
iso_to_cluster = {iso: cluster for cluster, members in cluster_map.items() for iso in members}

SHAPEFILE = "ne_110m_admin_0_countries.shp"

gdf_country = (
    gpd.read_file(SHAPEFILE)
    .loc[lambda d: d["ISO_A3"].isin(ISO38)]
    .assign(GeoCluster=lambda d: d["ISO_A3"].map(iso_to_cluster))
    .dropna(subset=["GeoCluster"])
)

country_cluster_lookup = gdf_country[["SOVEREIGNT", "ISO_A3", "GeoCluster"]].drop_duplicates()
gdf = gdf_country.dissolve(by="GeoCluster").reset_index()

regions = regionmask.Regions(
    outlines=gdf.geometry.values,
    names=gdf.GeoCluster.values,
    numbers=range(len(gdf))
)

mask2d = regions.mask(lon_or_obj=tg["longitude"], lat=tg["latitude"])

cluster_means = (
    xr.concat(
        [tg.where(mask2d == i).mean(dim=("latitude", "longitude")) for i in range(len(gdf))],
        dim="GeoCluster"
    )
    .assign_coords(GeoCluster=gdf.GeoCluster.values)
)

out_df = cluster_means.to_dataframe(name="tg_Celsius").reset_index()

out_with_country = country_cluster_lookup.merge(out_df, on="GeoCluster", how="left")

out_with_country.to_csv("GeoCluster_temp_with_countries.tsv", sep="\t", index=False)

df_daily = (
    out_with_country
    .rename(columns={"SOVEREIGNT": "Country", "tg_Celsius": "temperature"})
    .loc[:, ["Country", "GeoCluster", "time", "temperature"]]
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
    .agg(temperature=("temperature", "mean"))
)

avg_country = (
    monthly_country
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(temperature=("temperature", "mean"))
)

season_country = (
    monthly_country
    .groupby(["Country", "GeoCluster", "year"])
    .agg(temp_sd=("temperature", "std"))
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(temperature_seasonality=("temp_sd", "mean"))
)

avg_geo = (
    avg_country
    .groupby("GeoCluster", as_index=False)
    .agg(temperature=("temperature", "mean"))
)

print("Country-level temperature mean:")
print(avg_country.head(3))

print("\nCountry-level temperature seasonality:")
print(season_country.head(3))

print("\nGeoCluster-level temperature mean:")
print(avg_geo)

geo_out = avg_geo.copy()

geo_out[["GeoCluster", "temperature"]].to_csv(
    "GeoCluster_temperature.tsv",
    sep="\t",
    index=False
)

print(geo_out)
