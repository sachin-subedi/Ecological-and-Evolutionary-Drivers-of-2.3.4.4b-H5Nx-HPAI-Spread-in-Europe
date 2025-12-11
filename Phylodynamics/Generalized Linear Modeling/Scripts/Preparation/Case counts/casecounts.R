### Load libraries
library(dplyr)
library(readxl)
library(lubridate)
library(readr)
library(tibble)
library(tidyr)

  setwd("/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/GLM_code_predictors/codes/case_counts/")
case_data <- read_excel("Latest Reported Events.xlsx")

case_data <- case_data %>%
  mutate(Year = year(as.Date(`observation date`)))

case_data <- case_data %>%
  mutate(Country = recode(Country,
                          "U.K. of Great Britain and Northern Ireland" = "United Kingdom",
                          "Bosnia and Herzegovina" = "Bosnia Herz",
                          "Czech Republic" = "Czechia"
  ))

valid_countries <- c(
  "Netherlands", "Germany", "Czechia", "Belgium", "Denmark", "Sweden", "Moldova",
  "Slovenia", "Bulgaria", "Italy", "Ukraine", "United Kingdom", "Ireland",
  "Spain", "Croatia", "Poland", "Romania", "Albania", "Kosovo", "Greece", "France",
  "Norway", "Luxembourg", "Iceland", "Latvia", "Lithuania", "Finland", "Cyprus",
  "Austria", "Slovakia", "Switzerland", "Estonia", "Bosnia Herz", "Serbia",
  "Portugal", "North Macedonia", "Hungary"
)

case_data <- case_data %>%
  filter(Country %in% valid_countries)

cluster_region_lookup <- tribble(
  ~Country,                  ~GeoCluster,         ~SeqCluster,
  "Albania",                 "GeoCluster_One",     "SeqCluster_Three",
  "Austria",                 "GeoCluster_Five",    "SeqCluster_Fifteen",
  "Belgium",                 "GeoCluster_Three",   "SeqCluster_Seven",
  "Bosnia and Herzegovina",  "GeoCluster_Five",    "SeqCluster_Twelve",
  "Bulgaria",                "GeoCluster_One",     "SeqCluster_Fourteen",
  "Croatia",                 "GeoCluster_Five",    "SeqCluster_Twelve",
  "Cyprus",                  "GeoCluster_One",     "SeqCluster_Fourteen",
  "Czechia",                 "GeoCluster_Five",   "SeqCluster_Fifteen",
  "Denmark",                 "GeoCluster_Three",    "SeqCluster_Eleven",
  "Estonia",                 "GeoCluster_Four",     "SeqCluster_Four",
  "Finland",                 "GeoCluster_Four",     "SeqCluster_Four",
  "France",                  "GeoCluster_Three",     "SeqCluster_Two",
  "Germany",                 "GeoCluster_Three",     "SeqCluster_Nine",
  "Greece",                  "GeoCluster_One",    "SeqCluster_Three",
  "Hungary",                 "GeoCluster_Two",   "SeqCluster_Five",
  "Iceland",                 "GeoCluster_Three",    "SeqCluster_One",
  "Ireland",                 "GeoCluster_Three",     "SeqCluster_One",
  "Italy",                   "GeoCluster_Five",     "SeqCluster_Eight",
  "Kosovo",                  "GeoCluster_One",   "SeqCluster_Three",
  "Latvia",                  "GeoCluster_Four",     "SeqCluster_Four",
  "Lithuania",               "GeoCluster_Four",     "SeqCluster_Four",
  "Luxembourg",              "GeoCluster_Three",   "SeqCluster_Seven",
  "Moldova",                 "GeoCluster_Two",    "SeqCluster_Five",
  "Netherlands",             "GeoCluster_Three",   "SeqCluster_Seven",
  "North Macedonia",         "GeoCluster_One",    "SeqCluster_Three",
  "Norway",                  "GeoCluster_Four",    "SeqCluster_Thirteen",
  "Poland",                  "GeoCluster_Two",     "SeqCluster_Five",
  "Portugal",                "GeoCluster_Three",     "SeqCluster_Ten",
  "Romania",                 "GeoCluster_One",     "SeqCluster_Fourteen",
  "Serbia",                  "GeoCluster_One",    "SeqCluster_Three",
  "Slovakia",                "GeoCluster_Two",   "SeqCluster_Five",
  "Slovenia",                "GeoCluster_Five",    "SeqCluster_Fifteen",
  "Spain",                   "GeoCluster_Three",     "SeqCluster_Ten",
  "Sweden",                  "GeoCluster_Four",     "SeqCluster_SeqCluster_",
  "Switzerland",             "GeoCluster_Three",    "SeqCluster_Two",
  "Ukraine",                 "GeoCluster_Two",     "SeqCluster_Five",
  "United Kingdom",          "GeoCluster_Three",   "SeqCluster_Six"
)

case_data <- case_data %>%
  left_join(cluster_region_lookup, by = "Country")

GeoCluster_counts <- case_data %>%
  count(GeoCluster, name = "case_counts")

# Ensure all GeoClusters are represented
all_geo <- tibble(GeoCluster = unique(cluster_region_lookup$GeoCluster))
GeoCluster_counts <- all_geo %>%
  left_join(GeoCluster_counts, by = "GeoCluster") %>%
  mutate(case_counts = as.numeric(case_counts),
         case_counts = replace_na(case_counts, 0.000001),
         case_counts = format(case_counts, scientific = FALSE, trim = TRUE)) %>%
  arrange(desc(as.numeric(case_counts)))

write_tsv(GeoCluster_counts, "GeoCluster_case_counts.tsv")

