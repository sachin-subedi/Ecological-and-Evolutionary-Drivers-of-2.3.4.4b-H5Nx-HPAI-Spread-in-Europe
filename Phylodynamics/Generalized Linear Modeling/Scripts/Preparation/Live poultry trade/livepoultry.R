# ── 1. Load Required Libraries ────────────────────────────
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

setwd("livepoultry/")
input_file <- "TradeData_4_17_2025_0_38_32.csv"

iso3_map <- c(
  "NLD" = "Netherlands", "DEU" = "Germany", "CZE" = "Czechia",
  "BEL" = "Belgium", "DNK" = "Denmark", "SWE" = "Sweden", "MDA" = "Moldova",
  "SVN" = "Slovenia", "BGR" = "Bulgaria", "ITA" = "Italy", "UKR" = "Ukraine",
  "GBR" = "United Kingdom", "IRL" = "Ireland", "ESP" = "Spain", "HRV" = "Croatia",
  "POL" = "Poland", "ROU" = "Romania", "ALB" = "Albania", "XKX" = "Kosovo",
  "GRC" = "Greece", "FRA" = "France", "NOR" = "Norway", "LUX" = "Luxembourg",
  "ISL" = "Iceland", "LVA" = "Latvia", "LTU" = "Lithuania", "FIN" = "Finland",
  "CYP" = "Cyprus", "AUT" = "Austria", "SVK" = "Slovakia", "CHE" = "Switzerland",
  "EST" = "Estonia", "BIH" = "Bosnia and Herzegovina", "SRB" = "Republic of Serbia",
  "PRT" = "Portugal", "MKD" = "North Macedonia", "HUN" = "Hungary"
)

selected_iso <- names(iso3_map)

cluster_region_lookup <- tribble(
  ~Country,                  ~GeoCluster,
  "Albania","HC_Cluster_2_Mediterranean",
  "Austria","HC_Cluster_2_Alpine",
  "Belgium","HC_Cluster_1_Atlantic",
  "Bosnia and Herzegovina","HC_Cluster_2_Alpine",
  "Bulgaria","HC_Cluster_2_Continental",
  "Croatia","HC_Cluster_2_Continental",
  "Cyprus","HC_Cluster_2_Mediterranean",
  "Czechia","HC_Cluster_2_Continental",
  "Denmark","HC_Cluster_2_Continental",
  "Estonia","HC_Cluster_3_Boreal",
  "Finland","HC_Cluster_3_Boreal",
  "France","HC_Cluster_1_Atlantic",
  "Germany","HC_Cluster_1_Continental",
  "Greece","HC_Cluster_2_Mediterranean",
  "Hungary","HC_Cluster_2_Pannonian",
  "Iceland","HC_Cluster_1_Atlantic",
  "Ireland","HC_Cluster_1_Atlantic",
  "Italy","HC_Cluster_2_Mediterranean",
  "Kosovo","HC_Cluster_2_Continental",
  "Latvia","HC_Cluster_3_Boreal",
  "Lithuania","HC_Cluster_3_Boreal",
  "Luxembourg","HC_Cluster_1_Continental",
  "Moldova","HC_Cluster_2_Continental",
  "Netherlands","HC_Cluster_1_Atlantic",
  "North Macedonia","HC_Cluster_2_Continental",
  "Norway","HC_Cluster_3_Alpine",
  "Poland","HC_Cluster_2_Continental",
  "Portugal","HC_Cluster_4_Mediterranean",
  "Romania","HC_Cluster_2_Continental",
  "Republic of Serbia","HC_Cluster_2_Continental",
  "Slovakia","HC_Cluster_2_Alpine",
  "Slovenia","HC_Cluster_2_Continental",
  "Spain","HC_Cluster_4_Mediterranean",
  "Sweden","HC_Cluster_3_Boreal",
  "Switzerland","HC_Cluster_1_Alpine",
  "Ukraine","HC_Cluster_2_Continental",
  "United Kingdom","HC_Cluster_1_Atlantic"
)

DT <- fread(input_file, encoding = "Latin-1", header = TRUE, quote = "\"", fill = TRUE, showProgress = FALSE)

blank_cols <- names(DT)[vapply(DT, function(x) all(is.na(x)), logical(1))]
if (length(blank_cols)) DT[, (blank_cols) := NULL]

if ("row.names" %in% names(DT)) DT[, row.names := NULL]

# Convert types
DT[, Year := as.numeric(substr(refPeriodId, 1, 4))]
DT[, primaryValue := as.numeric(primaryValue)]
DT[, netWgt := as.numeric(netWgt)]

species_pattern <- "(fowl|gallus|duck|goose|geese|turkey|guinea)"

DT_filtered <- DT[
  Year >= 2016 & Year <= 2025 &
    grepl("live", cmdDesc, ignore.case = TRUE) &
    grepl(species_pattern, cmdDesc, ignore.case = TRUE) &
    reporterISO %in% selected_iso &
    partnerISO %in% selected_iso
]

DT_filtered <- DT_filtered %>%
  mutate(
    origin = iso3_map[reporterISO],
    destination = iso3_map[partnerISO]
  ) %>%
  left_join(cluster_region_lookup, by = c("origin" = "Country")) %>%
  rename(origin_GeoCluster = GeoCluster) %>%
  left_join(cluster_region_lookup, by = c("destination" = "Country")) %>%
  rename(dest_GeoCluster = GeoCluster)

missing_check <- DT_filtered %>%
  filter(is.na(origin_GeoCluster) | is.na(dest_GeoCluster)) %>%
  distinct(origin, destination)

if (nrow(missing_check) > 0) {
  print("⚠️ Missing cluster mappings:")
  print(missing_check)
}

make_matrix <- function(df, row_var, col_var, value_var) {
  df %>%
    filter(!is.na(.data[[row_var]]), !is.na(.data[[col_var]])) %>%
    group_by(across(all_of(c(row_var, col_var)))) %>%
    summarise(value = sum(.data[[value_var]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(col_var), values_from = value, values_fill = 0) %>%
    column_to_rownames(var = row_var) %>%
    {
      m <- as.matrix(.)
      diag(m) <- 0.000001
      m[m == 0] <- 0.000001
      as.data.frame(m)
    } %>%
    rownames_to_column(var = row_var)
}

geo_matrix_value  <- make_matrix(DT_filtered, "origin_GeoCluster", "dest_GeoCluster", "primaryValue")
geo_matrix_weight <- make_matrix(DT_filtered, "origin_GeoCluster", "dest_GeoCluster", "netWgt")

options(scipen = 999)

write_tsv(geo_matrix_value,  "live_poultry_trade_matrix_GeoCluster_value.tsv")
write_tsv(geo_matrix_weight, "live_poultry_trade_matrix_GeoCluster_weight.tsv")