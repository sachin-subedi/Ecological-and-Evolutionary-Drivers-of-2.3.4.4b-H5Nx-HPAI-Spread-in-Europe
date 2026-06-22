setwd("combined_traits/results/rates/Combined/")

library(dplyr)
library(readr)

region_labels <- c(
  "HC1_Alp" = "Central Alpine",
  "HC1_Atl" = "Atlantic",
  "HC1_Con" = "Western Continental",
  "HC2_Alp" = "Eastern Alpine",
  "HC2_Con" = "Eastern Continental",
  "HC2_Med" = "Southeast Mediterranean",
  "HC2_Pan" = "Pannonian",
  "HC3_Alp" = "Scandinavian Highlands",
  "HC3_Bor" = "Boreal Baltic",
  "HC4_Med" = "Iberian"
)

hab_map <- c(
  CM = "Coastal",
  FA = "Farm",
  FO = "Forest",
  GW = "Grassland",
  WT = "Wetland",
  UB = "Urban"
)

recode_label <- function(x) {
  x_chr <- as.character(x)
  sapply(x_chr, function(val) {
    parts    <- strsplit(val, "_")[[1]]
    hab_code <- tail(parts, 1)                        # e.g. "WT"
    geo_code <- paste(parts[-length(parts)], collapse = "_")  # e.g. "HC1_Atl"
    geo_full <- region_labels[geo_code]
    hab_full <- hab_map[hab_code]
    if (!is.na(geo_full) && !is.na(hab_full)) {
      paste(geo_full, hab_full)                       # "Atlantic Wetland"
    } else {
      val
    }
  }, USE.NAMES = FALSE)
}

all_levels <- as.vector(outer(unname(region_labels), unname(hab_map), paste))

df_equal <- read_csv("HG_bf_equal.csv") %>%
  mutate(subsample = "equal")

df_proportional <- read_csv("HG_bf_proportional.csv") %>%
  mutate(subsample = "proportional")

df_stratified <- read_csv("HG_bf_stratified.csv") %>%
  mutate(subsample = "stratified")

df_combined <- bind_rows(df_equal, df_proportional, df_stratified) %>%
  mutate(
    from      = recode_label(from),
    to        = recode_label(to),
    from      = factor(from, levels = all_levels),
    to        = factor(to,   levels = all_levels),
    subsample = factor(subsample, levels = c("equal", "proportional", "stratified"))
  )

write_csv(df_combined, "HG_bf_combined.csv")

cat("Done! Combined dimensions:", nrow(df_combined), "rows x", ncol(df_combined), "cols\n")
