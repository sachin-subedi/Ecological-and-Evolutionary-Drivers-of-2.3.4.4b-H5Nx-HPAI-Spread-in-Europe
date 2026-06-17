library(tidyverse)
library(sf)
library(rnaturalearth)

setwd("wildbird_counts")
csv_path <- "filtered_data_combined_Accipitridae.csv"

birds <- read_csv(
  csv_path,
  col_types = cols(
    startDayOfYear   = col_double(),
    individualCount  = col_double(),
    decimalLatitude  = col_double(),
    decimalLongitude = col_double(),
    family           = col_character(),
    genus            = col_character()
  )
) |>
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude),
         !is.na(individualCount)) |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs    = 4326,
           remove = FALSE)

europe <- rnaturalearth::ne_countries(continent = "Europe",
                                      returnclass = "sf") %>%
  filter(name != "Russia") %>%
  rename(Country = name)

birds <- st_join(
  birds,
  dplyr::select(europe, Country),
  join = st_within,
  left = FALSE
)

print(colnames(birds))
table(is.na(birds$Country))

cluster_region_lookup <- tribble(
  ~Country,                   ~GeoCluster,
  "Albania",                  "HC_Cluster_2_Mediterranean",
  "Austria",                  "HC_Cluster_2_Alpine",
  "Belarus",                  "HC_Cluster_2_Alpine",
  "Belgium",                  "HC_Cluster_1_Atlantic",
  "Bosnia and Herz.",         "HC_Cluster_2_Alpine",
  "Bulgaria",                 "HC_Cluster_2_Continental",
  "Croatia",                  "HC_Cluster_2_Continental",
  "Cyprus",                   "HC_Cluster_2_Mediterranean",
  "Czechia",                  "HC_Cluster_2_Continental",
  "Denmark",                  "HC_Cluster_2_Continental",
  "Estonia",                  "HC_Cluster_3_Boreal",
  "Finland",                  "HC_Cluster_3_Boreal",
  "France",                   "HC_Cluster_1_Atlantic",
  "Germany",                  "HC_Cluster_1_Continental",
  "Greece",                   "HC_Cluster_2_Mediterranean",
  "Hungary",                  "HC_Cluster_2_Pannonian",
  "Iceland",                  "HC_Cluster_1_Atlantic",
  "Ireland",                  "HC_Cluster_1_Atlantic",
  "Italy",                    "HC_Cluster_2_Mediterranean",
  "Kosovo",                   "HC_Cluster_2_Continental",
  "Latvia",                   "HC_Cluster_3_Boreal",
  "Lithuania",                "HC_Cluster_3_Boreal",
  "Luxembourg",               "HC_Cluster_1_Continental",
  "Moldova",                  "HC_Cluster_2_Continental",
  "Montenegro",          "HC_Cluster_2_Continental",
  "Norway",                   "HC_Cluster_2_Alpine",
  "Netherlands",              "HC_Cluster_1_Atlantic",
  "North Macedonia",          "HC_Cluster_2_Continental",
  "Norway",                   "HC_Cluster_3_Alpine",
  "Poland",                   "HC_Cluster_2_Continental",
  "Portugal",                 "HC_Cluster_4_Mediterranean",
  "Romania",                  "HC_Cluster_2_Continental",
  "Serbia",                   "HC_Cluster_2_Continental",
  "Republic of Serbia",       "HC_Cluster_2_Continental",
  "Slovakia",                 "HC_Cluster_2_Alpine",
  "Slovenia",                 "HC_Cluster_2_Continental",
  "Spain",                    "HC_Cluster_4_Mediterranean",
  "Sweden",                   "HC_Cluster_3_Boreal",
  "Switzerland",              "HC_Cluster_1_Alpine",
  "Ukraine",                  "HC_Cluster_2_Continental",
  "United Kingdom",           "HC_Cluster_1_Atlantic"
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

birds_gc <- birds %>%
  left_join(cluster_region_lookup, by = "Country") %>%
  filter(!is.na(GeoCluster))

unmatched <- birds %>%
  st_drop_geometry() %>%
  distinct(Country) %>%
  anti_join(cluster_region_lookup, by = "Country")
if (nrow(unmatched) > 0) message("⚠️  Unmatched countries: ", paste(unmatched$Country, collapse = ", "))

counts_year_cluster <- birds_gc %>%
  st_drop_geometry() %>%
  filter(year %in% 2016:2025) %>%
  group_by(GeoCluster, year) %>%
  summarise(total_count = sum(individualCount, na.rm = TRUE),
            .groups = "drop")

counts_complete <- counts_year_cluster %>%
  tidyr::complete(GeoCluster = hc_clusters,   # enforce all 10 clusters present
                  year = 2016:2025,
                  fill = list(total_count = 0)) %>%
  arrange(GeoCluster, year)

cluster_mean <- counts_complete %>%
  group_by(GeoCluster) %>%
  summarise(Accipitridae_counts = mean(total_count),
            .groups = "drop")

print(cluster_mean, n = Inf)

out_dir <- dirname(csv_path)

write_tsv(counts_complete,
          file.path(out_dir, "HC_Accipitridae_counts_by_GeoCluster_year.tsv"))

write_tsv(cluster_mean,
          file.path(out_dir, "HC_Accipitridae_Accipitridae_counts.tsv"))

