##Temperature
import xarray as xr
import rioxarray
import geopandas as gpd
import regionmask
import pandas as pd
import xarray as xr

ds = xr.open_dataset(
    "tg_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc",
    engine="netcdf4",      # 👈 force the NetCDF4 backend
    chunks={"time": 365}
)

ds = ds.rio.write_crs(4326)
tg = ds["tg"]
print(ds)

# ── 2. Define countries and GeoClusters ─────────────────────────────────────
ISO38 = [
    "ALB","AUT","BEL","BIH","BGR","CHE","CYP","CZE","DEU","DNK","ESP","EST",
    "FIN","FRA","GBR","GRC","HRV","HUN","IRL","ISL","ITA","XKX","LTU","LUX",
    "LVA","MDA","MKD","MNE","NLD","NOR","POL","PRT","ROU","SRB","SVK","SVN",
    "SWE","UKR","BLR"
]

cluster_map = {
    "One":   ["ALB", "BGR", "CYP", "GRC", "XKX", "MKD", "SRB", "ROU"],
    "Two":   ["HUN", "SVK", "POL", "UKR", "MDA", "BLR"],
    "Three": ["DNK", "FRA", "BEL", "DEU", "ISL", "IRL", "LUX", "NLD", "PRT", "ESP", "CHE", "GBR"],
    "Four":  ["EST", "FIN", "LVA", "LTU", "NOR", "SWE"],
    "Five":  ["AUT", "MNE", "BIH", "HRV", "CZE", "ITA", "SVN"]
}
iso_to_cluster = {iso: cluster for cluster, members in cluster_map.items() for iso in members}

# ── 3. Read shapefile and assign cluster and country ────────────────────────
SHAPEFILE = "/home/ss11645/Landcover/shapes/ne_110m_admin_0_countries.shp"

# Read full country-level geometry
gdf_country = (
    gpd.read_file(SHAPEFILE)
    .loc[lambda d: d["ISO_A3"].isin(ISO38)]
    .assign(GeoCluster=lambda d: d["ISO_A3"].map(iso_to_cluster))
    .dropna(subset=["GeoCluster"])
)

# Save country–cluster mapping with names
country_cluster_lookup = gdf_country[["SOVEREIGNT", "ISO_A3", "GeoCluster"]].drop_duplicates()
gdf = gdf_country.dissolve(by="GeoCluster").reset_index()

# ── 5. Build mask from merged cluster geometries ────────────────────────────
regions = regionmask.Regions(
    outlines=gdf.geometry.values,
    names=gdf.GeoCluster.values,
    numbers=range(len(gdf))
)

mask2d = regions.mask(lon_or_obj=tg["longitude"], lat=tg["latitude"])

# ── 6. Compute average temperature per cluster per day ──────────────────────
cluster_means = (
    xr.concat(
        [tg.where(mask2d == i).mean(dim=("latitude", "longitude")) for i in range(len(gdf))],
        dim="GeoCluster"
    )
    .assign_coords(GeoCluster=gdf.GeoCluster.values)
)

# ── 7. Convert to DataFrame ─────────────────────────────────────────────────
out_df = cluster_means.to_dataframe(name="tg_Celsius").reset_index()

# ── 8. Merge country names with GeoCluster average ──────────────────────────
out_with_country = country_cluster_lookup.merge(out_df, on="GeoCluster", how="left")

# ── 9. Export final result ──────────────────────────────────────────────────
out_with_country.to_csv("GeoCluster_temp_with_countries.tsv", sep="\t", index=False)

# ── 1. Rename and parse date ───────────────────────────────────────────────
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

# ── 2. DAILY → MONTHLY MEAN (per Country) ───────────────────────────────────
monthly_country = (
    df_daily
    .groupby(["Country", "GeoCluster", "year", "month"], as_index=False)
    .agg(temperature=("temperature", "mean"))
)

# ── 3. MULTI-YEAR COUNTRY MEAN ─────────────────────────────────────────────
avg_country = (
    monthly_country
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(temperature=("temperature", "mean"))
)

# ── 4. TEMPERATURE SEASONALITY PER COUNTRY ─────────────────────────────────
season_country = (
    monthly_country
    .groupby(["Country", "GeoCluster", "year"])
    .agg(temp_sd=("temperature", "std"))
    .groupby(["Country", "GeoCluster"], as_index=False)
    .agg(temperature_seasonality=("temp_sd", "mean"))
)

# ── 5. AGGREGATE TO GEOCLUSTER LEVEL ───────────────────────────────────────
avg_geo = (
    avg_country
    .groupby("GeoCluster", as_index=False)
    .agg(temperature=("temperature", "mean"))
)

# ── 7. View summaries in memory ────────────────────────────────────────────
print("Country-level temperature mean:")
print(avg_country.head(3))

print("\nCountry-level temperature seasonality:")
print(season_country.head(3))

print("\nGeoCluster-level temperature mean:")
print(avg_geo)

# ── 8. Add prefix and export GeoCluster-level temperature to TSV ────────────
geo_out = avg_geo.copy()
geo_out["GeoCluster"] = "GeoCluster_" + geo_out["GeoCluster"].astype(str)

geo_out[["GeoCluster", "temperature"]].to_csv(
    "GeoCluster_temperature.tsv",
    sep="\t",
    index=False
)

print(geo_out)
