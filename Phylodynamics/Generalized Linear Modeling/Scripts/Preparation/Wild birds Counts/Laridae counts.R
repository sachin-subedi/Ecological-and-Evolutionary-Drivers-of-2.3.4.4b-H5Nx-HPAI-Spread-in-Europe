library(tidyverse)
library(sf)
library(rnaturalearth)

csv_path <- "filtered_data_combined_Laridae.csv"

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
           crs = 4326,
           remove = FALSE)



library(dplyr)

europe <- rnaturalearth::ne_countries(continent = "Europe",
                                      returnclass = "sf") %>%
  filter(name != "Russia") %>%
  rename(Country = name)

# spatial join
birds <- st_join(
  birds,
  dplyr::select(europe, Country),
  join = st_within,
  left = FALSE
)

print(colnames(birds))

table(is.na(birds$Country))
cluster_region_lookup <- tribble(
  ~Country,                  ~GeoCluster,
  "Albania",                 "GeoCluster_One",
  "Austria",                 "GeoCluster_Five",
  "Belgium",                 "GeoCluster_Three",
  "Bosnia and Herzegovina",  "GeoCluster_Five",
  "Bulgaria",                "GeoCluster_One",
  "Croatia",                 "GeoCluster_Five",
  "Cyprus",                  "GeoCluster_One",
  "Czechia",                 "GeoCluster_Five",
  "Denmark",                 "GeoCluster_Three",
  "Estonia",                 "GeoCluster_Four",
  "Finland",                 "GeoCluster_Four",
  "France",                  "GeoCluster_Three",
  "Germany",                 "GeoCluster_Three",
  "Greece",                  "GeoCluster_One",
  "Hungary",                 "GeoCluster_Two",
  "Iceland",                 "GeoCluster_Three",
  "Ireland",                 "GeoCluster_Three",
  "Italy",                   "GeoCluster_Five",
  "Kosovo",                  "GeoCluster_One",
  "Latvia",                  "GeoCluster_Four",
  "Lithuania",               "GeoCluster_Four",
  "Luxembourg",              "GeoCluster_Three",
  "Moldova",                 "GeoCluster_Two",
  "Netherlands",             "GeoCluster_Three",
  "North Macedonia",         "GeoCluster_One",
  "Norway",                  "GeoCluster_Four",
  "Poland",                  "GeoCluster_Two",
  "Portugal",                "GeoCluster_Three",
  "Romania",                 "GeoCluster_One",
  "Serbia",                  "GeoCluster_One",
  "Slovakia",                "GeoCluster_Two",
  "Slovenia",                "GeoCluster_Five",
  "Spain",                   "GeoCluster_Three",
  "Sweden",                  "GeoCluster_Four",
  "Switzerland",             "GeoCluster_Three",
  "Ukraine",                 "GeoCluster_Two",
  "United Kingdom",          "GeoCluster_Three"
)

library(tidyverse)

birds_gc <- birds %>%
  left_join(cluster_region_lookup, by = "Country") %>%
  filter(!is.na(GeoCluster))

counts_year_cluster <- birds_gc %>%
  st_drop_geometry() %>%
  filter(year %in% 2016:2025) %>%
  group_by(GeoCluster, year) %>%
  summarise(total_count = sum(individualCount, na.rm = TRUE),
            .groups = "drop")

counts_complete <- counts_year_cluster %>%
  tidyr::complete(GeoCluster, year = 2016:2025,
                  fill = list(total_count = 0)) %>%
  arrange(GeoCluster, year)

cluster_mean <- counts_complete %>%
  group_by(GeoCluster) %>%
  summarise(laridae_counts = mean(total_count),
            .groups = "drop")

out_dir <- dirname(csv_path)
write_tsv(counts_complete,
          file.path(out_dir, "Laridae_counts_by_GeoCluster_year.tsv"))

write_tsv(cluster_mean,
          file.path(out_dir, "Laridae_laridae_counts.tsv"))
 

