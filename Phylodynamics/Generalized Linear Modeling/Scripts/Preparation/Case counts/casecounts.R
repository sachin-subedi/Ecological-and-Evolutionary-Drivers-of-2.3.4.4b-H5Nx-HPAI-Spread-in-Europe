library(dplyr)
library(readxl)
library(lubridate)
library(readr)
library(tibble)
library(tidyr)

setwd("casecounts")
case_data <- read_excel("Latest Reported Events.xlsx")

case_data <- case_data %>%
  mutate(Year = year(as.Date(`observation date`)))

case_data <- case_data %>%
  mutate(Country = recode(Country,
                          "U.K. of Great Britain and Northern Ireland" = "United Kingdom",
                          "Czech Republic"                             = "Czechia",
                          "Bosnia and Herzegovina"                     = "Bosnia and Herz."   # matches HC lookup
  ))

valid_countries <- c(
  "Albania", "Austria", "Belarus", "Belgium", "Bosnia and Herz.", "Bulgaria",
  "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", "Finland", "France",
  "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Kosovo",
  "Latvia", "Lithuania", "Luxembourg", "Moldova", "Montenegro", "Netherlands",
  "North Macedonia", "Norway", "Poland", "Portugal", "Romania",
  "Serbia", "Republic of Serbia", "Slovakia", "Slovenia", "Spain",
  "Sweden", "Switzerland", "Ukraine", "United Kingdom"
)

case_data <- case_data %>%
  filter(Country %in% valid_countries)

cluster_region_lookup <- tribble(
  ~Country,                ~GeoCluster,
  "Albania",               "HC_Cluster_2_Mediterranean",
  "Austria",               "HC_Cluster_2_Alpine",
  "Belarus",               "HC_Cluster_2_Alpine",
  "Belgium",               "HC_Cluster_1_Atlantic",
  "Bosnia and Herz.",      "HC_Cluster_2_Alpine",
  "Bulgaria",              "HC_Cluster_2_Continental",
  "Croatia",               "HC_Cluster_2_Continental",
  "Cyprus",                "HC_Cluster_2_Mediterranean",
  "Czechia",               "HC_Cluster_2_Continental",
  "Denmark",               "HC_Cluster_2_Continental",
  "Estonia",               "HC_Cluster_3_Boreal",
  "Finland",               "HC_Cluster_3_Boreal",
  "France",                "HC_Cluster_1_Atlantic",
  "Germany",               "HC_Cluster_1_Continental",
  "Greece",                "HC_Cluster_2_Mediterranean",
  "Hungary",               "HC_Cluster_2_Pannonian",
  "Iceland",               "HC_Cluster_1_Atlantic",
  "Ireland",               "HC_Cluster_1_Atlantic",
  "Italy",                 "HC_Cluster_2_Mediterranean",
  "Kosovo",                "HC_Cluster_2_Continental",
  "Latvia",                "HC_Cluster_3_Boreal",
  "Lithuania",             "HC_Cluster_3_Boreal",
  "Luxembourg",            "HC_Cluster_1_Continental",
  "Moldova",               "HC_Cluster_2_Continental",
  "Montenegro",            "HC_Cluster_2_Continental",
  "Netherlands",           "HC_Cluster_1_Atlantic",
  "North Macedonia",       "HC_Cluster_2_Continental",
  "Norway",                "HC_Cluster_3_Alpine",
  "Poland",                "HC_Cluster_2_Continental",
  "Portugal",              "HC_Cluster_4_Mediterranean",
  "Romania",               "HC_Cluster_2_Continental",
  "Serbia",                "HC_Cluster_2_Continental",
  "Republic of Serbia",    "HC_Cluster_2_Continental",
  "Slovakia",              "HC_Cluster_2_Alpine",
  "Slovenia",              "HC_Cluster_2_Continental",
  "Spain",                 "HC_Cluster_4_Mediterranean",
  "Sweden",                "HC_Cluster_3_Boreal",
  "Switzerland",           "HC_Cluster_1_Alpine",
  "Ukraine",               "HC_Cluster_2_Continental",
  "United Kingdom",        "HC_Cluster_1_Atlantic"
) %>% as.data.frame()

hc_clusters <- c(
  "HC_Cluster_1_Alpine",
  "HC_Cluster_1_Atlantic",
  "HC_Cluster_1_Continental",
  "HC_Cluster_2_Alpine",
  "HC_Cluster_2_Continental",
  "HC_Cluster_2_Mediterranean",
  "HC_Cluster_2_Pannonian",
  "HC_Cluster_3_Alpine",
  "HC_Cluster_3_Boreal",
  "HC_Cluster_4_Mediterranean"
)

case_data <- case_data %>%
  left_join(cluster_region_lookup, by = "Country")

unmatched <- case_data %>%
  filter(is.na(GeoCluster)) %>%
  distinct(Country)
if (nrow(unmatched) > 0) {
  warning("⚠️  Unmatched countries (no HC cluster): ",
          paste(unmatched$Country, collapse = ", "))
}

GeoCluster_counts <- case_data %>%
  count(GeoCluster, name = "case_counts")

GeoCluster_counts <- tibble(GeoCluster = hc_clusters) %>%
  left_join(GeoCluster_counts, by = "GeoCluster") %>%
  mutate(
    case_counts = as.numeric(case_counts),
    case_counts = replace_na(case_counts, 0.000001),
    case_counts = format(case_counts, scientific = FALSE, trim = TRUE)
  ) %>%
  arrange(desc(as.numeric(case_counts)))

write_tsv(GeoCluster_counts, "HC_GeoCluster_case_counts.tsv")

summary_stats <- function(df, colname) {
  df %>%
    mutate(case_counts = as.numeric(case_counts)) %>%
    summarise(
      mean  = round(mean(case_counts),  2),
      sd    = round(sd(case_counts),    2),
      min   = min(case_counts),
      max   = max(case_counts),
      range = max - min
    ) %>%
    mutate(Type = colname)
}

GeoCluster_summary <- summary_stats(GeoCluster_counts, "HC_GeoCluster")
print(GeoCluster_summary)

