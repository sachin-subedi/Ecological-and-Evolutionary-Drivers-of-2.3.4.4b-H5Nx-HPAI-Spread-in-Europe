library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

setwd("/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/GLM_code_predictors/codes/livepoultry_trade/")
input_file <- "TradeData_4_17_2025_0_38_32.csv"

iso3_map <- c(
  "NLD" = "Netherlands",
  "DEU" = "Germany",
  "CZE" = "Czechia",
  "BEL" = "Belgium",
  "DNK" = "Denmark",
  "SWE" = "Sweden",
  "MDA" = "Moldova",
  "SVN" = "Slovenia",
  "BGR" = "Bulgaria",
  "ITA" = "Italy",
  "UKR" = "Ukraine",
  "GBR" = "United_Kingdom",
  "IRL" = "Ireland",
  "ESP" = "Spain",
  "HRV" = "Croatia",
  "POL" = "Poland",
  "ROU" = "Romania",
  "ALB" = "Albania",
  "XKX" = "Kosovo",
  "GRC" = "Greece",
  "FRA" = "France",
  "NOR" = "Norway",
  "LUX" = "Luxembourg",
  "ISL" = "Iceland",
  "LVA" = "Latvia",
  "LTU" = "Lithuania",
  "FIN" = "Finland",
  "CYP" = "Cyprus",
  "AUT" = "Austria",
  "SVK" = "Slovakia",
  "CHE" = "Switzerland",
  "EST" = "Estonia",
  "BIH" = "Bosnia_Herz",
  "SRB" = "Serbia",
  "PRT" = "Portugal",
  "MKD" = "North_Macedonia",
  "HUN" = "Hungary"
)
country_list <- unname(iso3_map)

cluster_region_lookup <- tribble(
  ~Country, ~GeoCluster,
  "Albania", "GeoCluster_One",
  "Austria", "GeoCluster_Five",
  "Belgium", "GeoCluster_Three",
  "Bosnia_Herz", "GeoCluster_Five",
  "Bulgaria", "GeoCluster_One",
  "Croatia", "GeoCluster_Five",
  "Cyprus", "GeoCluster_One",
  "Czechia", "GeoCluster_Five",
  "Denmark", "GeoCluster_Three",
  "Estonia", "GeoCluster_Four",
  "Finland", "GeoCluster_Four",
  "France", "GeoCluster_Three",
  "Germany", "GeoCluster_Three",
  "Greece", "GeoCluster_One",
  "Hungary", "GeoCluster_Two",
  "Iceland", "GeoCluster_Three",
  "Ireland", "GeoCluster_Three",
  "Italy", "GeoCluster_Five",
  "Kosovo", "GeoCluster_One",
  "Latvia", "GeoCluster_Four",
  "Lithuania", "GeoCluster_Four",
  "Luxembourg", "GeoCluster_Three",
  "Moldova", "GeoCluster_Two",
  "Netherlands", "GeoCluster_Three",
  "North_Macedonia", "GeoCluster_One",
  "Norway", "GeoCluster_Four",
  "Poland", "GeoCluster_Two",
  "Portugal", "GeoCluster_Three",
  "Romania", "GeoCluster_One",
  "Serbia", "GeoCluster_One",
  "Slovakia", "GeoCluster_Two",
  "Slovenia", "GeoCluster_Five",
  "Spain", "GeoCluster_Three",
  "Sweden", "GeoCluster_Four",
  "Switzerland", "GeoCluster_Three",
  "Ukraine", "GeoCluster_Two",
  "United_Kingdom", "GeoCluster_Three"
)

DT <- fread(input_file, encoding = "Latin-1", header = TRUE, quote = "\"", fill = TRUE, showProgress = FALSE)
blank_cols <- names(DT)[vapply(DT, function(x) all(is.na(x)), logical(1))]
if (length(blank_cols)) DT[, (blank_cols) := NULL]
if ("row.names" %in% names(DT)) DT[, row.names := NULL]
DT[, Year := as.numeric(substr(refPeriodId, 1, 4))]
DT[, primaryValue := as.numeric(primaryValue)]
DT[, netWgt := as.numeric(netWgt)]

species_pattern <- "(fowl|gallus|duck|goose|geese|turkey|guinea)"
selected_iso <- names(iso3_map)
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

make_matrix <- function(df, row_var, col_var, value_var) {
  df %>%
    group_by(across(all_of(c(row_var, col_var)))) %>%
    summarise(
      value = sum(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = all_of(col_var),
      values_from = value,
      values_fill = 0
    ) %>%
    column_to_rownames(var = row_var) %>%
    {
      m <- as.matrix(.); diag(m) <- 0.000001; m[m == 0] <- 0.000001; as.data.frame(m)
    } %>%
    tibble::rownames_to_column(var = row_var)
}

land_matrix_value  <- make_matrix(DT_filtered, "origin_GeoCluster", "dest_GeoCluster",  "primaryValue")
land_matrix_weight <- make_matrix(DT_filtered, "origin_GeoCluster", "dest_GeoCluster",  "netWgt")

options(scipen = 999)
write_tsv(land_matrix_value,  "live_poultry_trade_matrix_GeoCluster_value.tsv")
write_tsv(land_matrix_weight, "live_poultry_trade_matrix_GeoCluster_weight.tsv")

cat("✅ GeoCluster matrices exported.\n")
