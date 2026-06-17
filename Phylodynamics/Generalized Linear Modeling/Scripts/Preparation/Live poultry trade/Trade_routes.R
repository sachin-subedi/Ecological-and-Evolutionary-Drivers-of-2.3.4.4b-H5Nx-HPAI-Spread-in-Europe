library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(grid)

setwd("livepoultry/")
options(scipen = 999)

state_labels <- c(
  "HC_Cluster_1_Alpine"        = "HC1 Alpine",
  "HC_Cluster_1_Atlantic"      = "HC1 Atlantic",
  "HC_Cluster_1_Continental"   = "HC1 Continental",
  "HC_Cluster_2_Alpine"        = "HC2 Alpine",
  "HC_Cluster_2_Continental"   = "HC2 Continental",
  "HC_Cluster_2_Mediterranean" = "HC2 Mediterranean",
  "HC_Cluster_2_Pannonian"     = "HC2 Pannonian",
  "HC_Cluster_3_Alpine"        = "HC3 Alpine",
  "HC_Cluster_3_Boreal"        = "HC3 Boreal",
  "HC_Cluster_4_Mediterranean" = "HC4 Mediterranean"
)

region_pal <- c(
  "HC_Cluster_1_Alpine"        = "#0072B2",
  "HC_Cluster_1_Atlantic"      = "#CD5C5C",
  "HC_Cluster_1_Continental"   = "#BCBD22",
  "HC_Cluster_2_Alpine"        = "#26A69A",
  "HC_Cluster_2_Continental"   = "#CC79A7",
  "HC_Cluster_2_Mediterranean" = "lightslategray",
  "HC_Cluster_2_Pannonian"     = "#6B4C1B",
  "HC_Cluster_3_Alpine"        = "#7570B3",
  "HC_Cluster_3_Boreal"        = "#006D2C",
  "HC_Cluster_4_Mediterranean" = "#56B4E9"
)

iso3_map <- c(
  "NLD" = "Netherlands",   "DEU" = "Germany",         "CZE" = "Czechia",
  "BEL" = "Belgium",       "DNK" = "Denmark",          "SWE" = "Sweden",
  "MDA" = "Moldova",       "SVN" = "Slovenia",          "BGR" = "Bulgaria",
  "ITA" = "Italy",         "UKR" = "Ukraine",           "GBR" = "United Kingdom",
  "IRL" = "Ireland",       "ESP" = "Spain",             "HRV" = "Croatia",
  "POL" = "Poland",        "ROU" = "Romania",           "ALB" = "Albania",
  "XKX" = "Kosovo",        "GRC" = "Greece",            "FRA" = "France",
  "NOR" = "Norway",        "LUX" = "Luxembourg",        "ISL" = "Iceland",
  "LVA" = "Latvia",        "LTU" = "Lithuania",         "FIN" = "Finland",
  "CYP" = "Cyprus",        "AUT" = "Austria",           "SVK" = "Slovakia",
  "CHE" = "Switzerland",   "EST" = "Estonia",           "BIH" = "Bosnia and Herzegovina",
  "SRB" = "Republic of Serbia", "PRT" = "Portugal",    "MKD" = "North Macedonia",
  "HUN" = "Hungary",       "BLR" = "Belarus"
)
selected_iso <- names(iso3_map)

cluster_region_lookup <- tribble(
  ~Country,                       ~GeoCluster,
  "Albania",                      "HC_Cluster_2_Mediterranean",
  "Austria",                      "HC_Cluster_2_Alpine",
  "Belarus",                      "HC_Cluster_2_Continental",
  "Belgium",                      "HC_Cluster_1_Atlantic",
  "Bosnia and Herzegovina",       "HC_Cluster_2_Alpine",
  "Bulgaria",                     "HC_Cluster_2_Continental",
  "Croatia",                      "HC_Cluster_2_Continental",
  "Cyprus",                       "HC_Cluster_2_Mediterranean",
  "Czechia",                      "HC_Cluster_2_Continental",
  "Denmark",                      "HC_Cluster_2_Continental",
  "Estonia",                      "HC_Cluster_3_Boreal",
  "Finland",                      "HC_Cluster_3_Boreal",
  "France",                       "HC_Cluster_1_Atlantic",
  "Germany",                      "HC_Cluster_1_Continental",
  "Greece",                       "HC_Cluster_2_Mediterranean",
  "Hungary",                      "HC_Cluster_2_Pannonian",
  "Iceland",                      "HC_Cluster_1_Atlantic",
  "Ireland",                      "HC_Cluster_1_Atlantic",
  "Italy",                        "HC_Cluster_2_Mediterranean",
  "Kosovo",                       "HC_Cluster_2_Continental",
  "Latvia",                       "HC_Cluster_3_Boreal",
  "Lithuania",                    "HC_Cluster_3_Boreal",
  "Luxembourg",                   "HC_Cluster_1_Continental",
  "Moldova",                      "HC_Cluster_2_Continental",
  "Montenegro",                     "HC_Cluster_2_Alpine",
  "Netherlands",                  "HC_Cluster_1_Atlantic",
  "North Macedonia",              "HC_Cluster_2_Continental",
  "Norway",                       "HC_Cluster_3_Alpine",
  "Poland",                       "HC_Cluster_2_Continental",
  "Portugal",                     "HC_Cluster_4_Mediterranean",
  "Romania",                      "HC_Cluster_2_Continental",
  "Republic of Serbia",           "HC_Cluster_2_Continental",
  "Slovakia",                     "HC_Cluster_2_Alpine",
  "Slovenia",                     "HC_Cluster_2_Continental",
  "Spain",                        "HC_Cluster_4_Mediterranean",
  "Sweden",                       "HC_Cluster_3_Boreal",
  "Switzerland",                  "HC_Cluster_1_Alpine",
  "Ukraine",                      "HC_Cluster_2_Continental",
  "United Kingdom",               "HC_Cluster_1_Atlantic"
)

message("• Loading trade data...")
DT <- fread(
  "TradeData_4_17_2025_0_38_32.csv",
  encoding = "Latin-1", header = TRUE,
  quote = "\"", fill = TRUE, showProgress = FALSE
)

blank_cols <- names(DT)[vapply(DT, function(x) all(is.na(x)), logical(1))]
if (length(blank_cols)) DT[, (blank_cols) := NULL]
if ("row.names" %in% names(DT)) DT[, row.names := NULL]

DT[, Year         := as.numeric(substr(refPeriodId, 1, 4))]
DT[, primaryValue := as.numeric(primaryValue)]
DT[, netWgt       := as.numeric(netWgt)]

species_pattern <- "(fowl|gallus|duck|goose|geese|turkey|guinea)"

DT_filtered <- DT[
  Year >= 2016 & Year <= 2025 &
    grepl("live",            cmdDesc, ignore.case = TRUE) &
    grepl(species_pattern,   cmdDesc, ignore.case = TRUE) &
    reporterISO %in% selected_iso &
    partnerISO  %in% selected_iso
] %>%
  mutate(
    origin      = iso3_map[reporterISO],
    destination = iso3_map[partnerISO]
  ) %>%
  left_join(cluster_region_lookup, by = c("origin"      = "Country")) %>%
  rename(origin_GeoCluster = GeoCluster) %>%
  left_join(cluster_region_lookup, by = c("destination" = "Country")) %>%
  rename(dest_GeoCluster = GeoCluster)

europe_sf <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(!name %in% c("Russia")) %>%
  rename(Country = name)

name_fix <- c(
  "Bosnia and Herz."  = "Bosnia and Herzegovina",
  "Serbia"            = "Republic of Serbia",
  "Macedonia"         = "North Macedonia",
  "Czech Rep."        = "Czechia",
  "Czechia"           = "Czechia"
)
europe_sf <- europe_sf %>%
  mutate(Country = recode(Country, !!!name_fix))

europe_data <- europe_sf %>%
  left_join(cluster_region_lookup, by = "Country") %>%
  mutate(GeoCluster = factor(GeoCluster, levels = names(state_labels)))

centroids <- europe_sf %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(Country, lon, lat)

edges <- DT_filtered %>%
  filter(!is.na(primaryValue), origin != destination) %>%
  group_by(origin, destination,
           origin_GeoCluster, dest_GeoCluster) %>%
  summarise(
    total_value  = sum(primaryValue, na.rm = TRUE),
    total_weight = sum(netWgt,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(centroids, by = c("origin"      = "Country")) %>%
  rename(lon_from = lon, lat_from = lat) %>%
  left_join(centroids, by = c("destination" = "Country")) %>%
  rename(lon_to = lon, lat_to = lat) %>%
  filter(!is.na(lon_from), !is.na(lat_from),
         !is.na(lon_to),   !is.na(lat_to))

threshold <- quantile(edges$total_value, 0.85, na.rm = TRUE)
edges_top <- edges %>%
  filter(total_value >= threshold) %>%
  mutate(flow_cat = case_when(
    total_value < quantile(total_value, 1/3) ~ "Low",
    total_value < quantile(total_value, 2/3) ~ "Medium",
    TRUE                                      ~ "High"
  )) %>%
  mutate(flow_cat = factor(flow_cat, levels = c("Low", "Medium", "High")))

p <- ggplot() +
  
  geom_sf(
    data      = europe_data,
    aes(fill  = GeoCluster),
    colour    = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values   = region_pal,
    labels   = state_labels,
    name     = "GeoCluster",
    na.value = "grey88",
    drop     = FALSE
  ) +
  
  geom_curve(
    data = edges_top,
    aes(
      x         = lon_from, y    = lat_from,
      xend      = lon_to,   yend = lat_to,
      linewidth = flow_cat
    ),
    colour    = "black",
    curvature = 0.25,
    alpha     = 0.65,
    arrow     = arrow(
      angle  = 18,
      type   = "open",
      length = unit(0.10, "inches")
    ),
    lineend = "round"
  ) +
  scale_linewidth_manual(
    name   = "Trade flow\n(top 15% by value)",
    values = c("Low" = 0.4, "Medium" = 1.0, "High" = 2.0)
  ) +
  
  coord_sf(xlim = c(-12, 45), ylim = c(30, 72), expand = FALSE) +
  
  labs(
    title    = "Live Poultry Trade Routes in Europe",
    subtitle = ""
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title        = element_text(face = "bold", size = 15),
    plot.subtitle     = element_text(size = 12, color = "grey40"),
    
    legend.position   = "right",
    legend.box        = "vertical",
    legend.box.just   = "left",
    legend.margin     = margin(0, 0, 0, 6),
    legend.spacing.y  = unit(0.3, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title      = element_text(face = "bold", size = 12),
    legend.text       = element_text(size = 11),
    legend.key.size   = unit(0.55, "cm"),
    
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.line   = element_blank(),
    panel.grid  = element_blank(),
    panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.5),

    plot.margin = margin(t = 8, r = 5, b = 8, l = 8)
  ) +
  guides(
    fill      = guide_legend(order = 1,
                             override.aes = list(colour = "white"),
                             keywidth = unit(0.6, "cm"),
                             keyheight = unit(0.55, "cm")),
    linewidth = guide_legend(order = 2,
                             override.aes = list(colour = "black"),
                             keywidth = unit(1.2, "cm"),
                             keyheight = unit(0.55, "cm"))
  )
p
ggsave(
  "LivePoultry_TradeRoutes_GeoCluster_v2.png",
  plot   = p,
  width  = 16,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

