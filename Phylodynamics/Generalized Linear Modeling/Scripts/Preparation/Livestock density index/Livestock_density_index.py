import pandas as pd

file_path = "estat_tai09.tsv"
df = pd.read_csv(file_path, sep='\t')

df[['freq', 'leg_form', 'so_eur', 'statinfo','unit', 'farmtype', 'geo_code']] = df.iloc[:, 0].str.split(',', expand=True)

df.columns = df.columns.str.strip()

year_cols = [col.strip() for col in df.columns if col.strip().isdigit() and 2016 <= int(col.strip()) <= 2024]
for col in year_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df['avg_2016_2024'] = df[year_cols].mean(axis=1, skipna=True)

df['country_code'] = df['geo_code'].str[:2]

country_sum = df.groupby('country_code', as_index=False)['avg_2016_2024'].sum()
country_sum.rename(columns={'avg_2016_2024': 'livestock_density_index'}, inplace=True)

country_sum = country_sum[~country_sum['country_code'].isin(['TR', 'MT'])]

country_lookup = pd.DataFrame([
    ("AT", "Austria"), ("BE", "Belgium"), ("BG", "Bulgaria"), ("CH", "Switzerland"),
    ("CY", "Cyprus"), ("CZ", "Czechia"), ("DE", "Germany"), ("DK", "Denmark"),
    ("EE", "Estonia"), ("EL", "Greece"), ("ES", "Spain"), ("FI", "Finland"),
    ("FR", "France"), ("HR", "Croatia"), ("HU", "Hungary"), ("IE", "Ireland"),
    ("IS", "Iceland"), ("IT", "Italy"), ("LT", "Lithuania"), ("LU", "Luxembourg"),
    ("LV", "Latvia"), ("ME", "Montenegro"), ("MK", "North Macedonia"),
    ("NL", "Netherlands"), ("PL", "Poland"), ("PT", "Portugal"), ("RO", "Romania"),
    ("RS", "Serbia"), ("SE", "Sweden"), ("SI", "Slovenia"), ("SK", "Slovakia"),
    ("UK", "United Kingdom"), ("XK", "Kosovo")
], columns=["country_code", "Country"])

merged_df = country_sum.merge(country_lookup, on="country_code", how="left")

cluster_map = {
    "GeoCluster_One": [
        "Albania", "Bulgaria", "Cyprus", "Greece", "Kosovo", "North Macedonia", "Serbia", "Romania"
    ],
    "GeoCluster_Two": [
        "Hungary", "Slovakia", "Poland", "Ukraine", "Moldova", "Belarus"
    ],
    "GeoCluster_Three": [
        "Denmark", "France", "Belgium", "Germany", "Iceland", "Ireland", "Luxembourg",
        "Netherlands", "Portugal", "Spain", "Switzerland", "United Kingdom"
    ],
    "GeoCluster_Four": [
        "Estonia", "Finland", "Latvia", "Lithuania", "Norway", "Sweden"
    ],
    "GeoCluster_Five": [
        "Austria", "Montenegro", "Bosnia and Herz.", "Croatia", "Czechia", "Italy", "Slovenia"
    ]
}

geo_cluster_lookup = []
for cluster, countries in cluster_map.items():
    for country in countries:
        geo_cluster_lookup.append((country, cluster))
geo_cluster_df = pd.DataFrame(geo_cluster_lookup, columns=["Country", "GeoCluster"])

merged_with_cluster = merged_df.merge(geo_cluster_df, on="Country", how="left")

cluster_summary = merged_with_cluster.groupby('GeoCluster', as_index=False)['livestock_density_index'].sum()
cluster_summary.rename(columns={'livestock_density_index': 'livestock_density_index'}, inplace=True)

cluster_summary.to_csv("GeoCluster_livestock_density_index.tsv", sep='\t', index=False)
