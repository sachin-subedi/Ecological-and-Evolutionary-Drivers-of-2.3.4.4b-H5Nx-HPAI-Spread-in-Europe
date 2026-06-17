import pandas as pd

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

iso3_to_iso2 = {
    "ALB": "AL", "AUT": "AT", "BEL": "BE", "BIH": "BA", "BGR": "BG",
    "CHE": "CH", "CYP": "CY", "CZE": "CZ", "DEU": "DE", "DNK": "DK",
    "ESP": "ES", "EST": "EE", "FIN": "FI", "FRA": "FR", "GBR": "UK",
    "GRC": "EL", "HRV": "HR", "HUN": "HU", "IRL": "IE", "ISL": "IS",
    "ITA": "IT", "XKX": "XK", "LTU": "LT", "LUX": "LU", "LVA": "LV",
    "MDA": "MD", "MKD": "MK", "MNE": "ME", "NLD": "NL", "NOR": "NO",
    "POL": "PL", "PRT": "PT", "ROU": "RO", "SRB": "RS", "SVK": "SK",
    "SVN": "SI", "SWE": "SE", "UKR": "UA", "BLR": "BY",
}

hc_cluster_lookup = [
    (iso3_to_iso2[iso3], cluster)
    for cluster, iso3_list in hc_cluster_map.items()
    for iso3 in iso3_list
    if iso3 in iso3_to_iso2
]
hc_cluster_df = pd.DataFrame(hc_cluster_lookup, columns=["country_code", "HC_Cluster"])

country_lookup = pd.DataFrame([
    ("AL", "Albania"),          ("AT", "Austria"),          ("BE", "Belgium"),
    ("BA", "Bosnia and Herzegovina"), ("BG", "Bulgaria"),   ("BY", "Belarus"),
    ("CH", "Switzerland"),      ("CY", "Cyprus"),            ("CZ", "Czechia"),
    ("DE", "Germany"),          ("DK", "Denmark"),           ("EE", "Estonia"),
    ("EL", "Greece"),           ("ES", "Spain"),             ("FI", "Finland"),
    ("FR", "France"),           ("HR", "Croatia"),           ("HU", "Hungary"),
    ("IE", "Ireland"),          ("IS", "Iceland"),           ("IT", "Italy"),
    ("LT", "Lithuania"),        ("LU", "Luxembourg"),        ("LV", "Latvia"),
    ("MD", "Moldova"),          ("ME", "Montenegro"),        ("MK", "North Macedonia"),
    ("NL", "Netherlands"),      ("NO", "Norway"),
    ("PL", "Poland"),           ("PT", "Portugal"),          ("RO", "Romania"),
    ("RS", "Serbia"),           ("SE", "Sweden"),            ("SI", "Slovenia"),
    ("SK", "Slovakia"),         ("UK", "United Kingdom"),    ("UA", "Ukraine"),
    ("XK", "Kosovo"),
], columns=["country_code", "Country"])

file_path = "estat_tgs00045.tsv"
df = pd.read_csv(file_path, sep='\t')

df[['freq', 'animals', 'unit', 'geo_code']] = df.iloc[:, 0].str.split(',', expand=True)

df.columns = df.columns.str.strip()

year_cols = [col.strip() for col in df.columns if col.strip().isdigit() and 2016 <= int(col.strip()) <= 2024]
for col in year_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df['avg_2016_2024'] = df[year_cols].mean(axis=1, skipna=True)

df['country_code'] = df['geo_code'].str.strip().str[:2]

country_sum = df.groupby('country_code', as_index=False)['avg_2016_2024'].sum()
country_sum.rename(columns={'avg_2016_2024': 'total_density'}, inplace=True)

country_sum = country_sum[~country_sum['country_code'].isin(['TR', 'MT'])]

if "NO" not in country_sum["country_code"].values:
    norway_row = pd.DataFrame([{"country_code": "NO", "total_density": 0.0000001}])
    country_sum = pd.concat([country_sum, norway_row], ignore_index=True)
    print("ℹ️  Norway ('NO') not found in Eurostat data — injected as 0.0000001")

country_sum = country_sum.merge(country_lookup, on="country_code", how="left")

merged = country_sum.merge(hc_cluster_df, on="country_code", how="left")

unmatched = merged[merged["HC_Cluster"].isna()]["country_code"].tolist()
if unmatched:
    print(f"⚠️  No HC cluster assigned for country_codes: {unmatched}")

hc_summary = (
    merged
    .dropna(subset=["HC_Cluster"])
    .groupby("HC_Cluster", as_index=False)["total_density"]
    .sum()
    .rename(columns={"total_density": "Animal_Population"})
)

output_file = "HC_GeoCluster_Animal_Population.tsv"
hc_summary.to_csv(output_file, sep='\t', index=False)
