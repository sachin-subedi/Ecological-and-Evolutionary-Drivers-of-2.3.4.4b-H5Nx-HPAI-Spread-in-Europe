library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(grid)
library(scales)
library(readr)

setwd("Region/")

options(scipen = 999)

state_labels <- c(
  "HC_Cluster_1_Alpine"        = "Central Alpine",
  "HC_Cluster_1_Atlantic"      = "Atlantic",
  "HC_Cluster_1_Continental"   = "Western Continental",
  "HC_Cluster_2_Alpine"        = "Eastern Alpine",
  "HC_Cluster_2_Continental"   = "Eastern Continental",
  "HC_Cluster_2_Mediterranean" = "Southeast Mediterranean",
  "HC_Cluster_2_Pannonian"     = "Pannonian",
  "HC_Cluster_3_Alpine"        = "Scandinavian Highlands",
  "HC_Cluster_3_Boreal"        = "Boreal Baltic",
  "HC_Cluster_4_Mediterranean" = "Iberian"
)

src_cols <- list(
  "HC_Cluster_1_Alpine"        = c(scales::alpha("#CC79A7",        0.25), "#CC79A7"),
  "HC_Cluster_1_Atlantic"      = c(scales::alpha("#CD5C5C",        0.25), "#CD5C5C"),
  "HC_Cluster_1_Continental"   = c(scales::alpha("#BCBD22",        0.25), "#BCBD22"),
  "HC_Cluster_2_Alpine"        = c(scales::alpha("#26A69A",        0.25), "#26A69A"),
  "HC_Cluster_2_Continental"   = c(scales::alpha("#0072B2",        0.25), "#0072B2"),
  "HC_Cluster_2_Mediterranean" = c(scales::alpha("lightslategray", 0.25), "lightslategray"),
  "HC_Cluster_2_Pannonian"     = c(scales::alpha("#6B4C1B",        0.25), "#6B4C1B"),
  "HC_Cluster_3_Alpine"        = c(scales::alpha("#7570B3",        0.25), "#7570B3"),
  "HC_Cluster_3_Boreal"        = c(scales::alpha("#006D2C",        0.25), "#006D2C"),
  "HC_Cluster_4_Mediterranean" = c(scales::alpha("#56B4E9",        0.25), "#56B4E9")
)

region_pal <- sapply(src_cols, `[[`, 2)

cluster_lookup <- data.frame(
  Country = c(
    "Albania", "Austria", "Belgium", "Bosnia and Herz.", "Bulgaria",
    "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia",
    "Finland", "France", "Germany", "Greece", "Hungary",
    "Ireland", "Italy", "Kosovo", "Latvia", "Lithuania",
    "Luxembourg", "Moldova", "Netherlands", "North Macedonia", "Norway",
    "Poland", "Portugal", "Romania", "Serbia", "Slovakia",
    "Slovenia", "Spain", "Sweden", "Switzerland", "Ukraine",
    "United Kingdom", "Belarus", "Montenegro"
  ),
  Region = c(
    "HC_Cluster_2_Mediterranean", "HC_Cluster_2_Alpine",      "HC_Cluster_1_Atlantic",
    "HC_Cluster_2_Alpine",        "HC_Cluster_2_Continental", "HC_Cluster_2_Continental",
    "HC_Cluster_2_Mediterranean", "HC_Cluster_2_Continental", "HC_Cluster_2_Continental",
    "HC_Cluster_3_Boreal",        "HC_Cluster_3_Boreal",      "HC_Cluster_1_Atlantic",
    "HC_Cluster_1_Continental",   "HC_Cluster_2_Mediterranean","HC_Cluster_2_Pannonian",
    "HC_Cluster_1_Atlantic",      "HC_Cluster_2_Mediterranean","HC_Cluster_2_Continental",
    "HC_Cluster_3_Boreal",        "HC_Cluster_3_Boreal",      "HC_Cluster_1_Continental",
    "HC_Cluster_2_Continental",   "HC_Cluster_1_Atlantic",    "HC_Cluster_2_Continental",
    "HC_Cluster_3_Alpine",        "HC_Cluster_2_Continental", "HC_Cluster_4_Mediterranean",
    "HC_Cluster_2_Continental",   "HC_Cluster_2_Continental", "HC_Cluster_2_Alpine",
    "HC_Cluster_2_Continental",   "HC_Cluster_4_Mediterranean","HC_Cluster_3_Boreal",
    "HC_Cluster_1_Alpine",        "HC_Cluster_2_Continental", "HC_Cluster_1_Atlantic",
    "HC_Cluster_2_Alpine",        "HC_Cluster_2_Alpine"
  ),
  stringsAsFactors = FALSE
)

europe <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(!name %in% c("Russia", "Iceland"))

europe_data <- europe %>%
  left_join(cluster_lookup, by = c("name" = "Country")) %>%
  mutate(Region = factor(Region, levels = names(state_labels)))

cluster_centers <- data.frame(
  Region = c(
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
  ),
  x = c(
    8.2,  # Central Alpine      — Switzerland
    -3.5, # Atlantic            — UK/France/Ireland/Belgium/NL
    10.5, # Western Continental — Germany/Luxembourg
    19.0, # Eastern Alpine      — Austria/Slovakia/Bosnia/Belarus/Montenegro
    28.0, # Eastern Continental — Poland/Romania/Ukraine/Serbia etc.
    22.0, # Southeast Mediterranean — Albania/Greece/Italy/Cyprus
    18.5, # Pannonian           — Hungary
    10.0, # Scandinavian Highlands — Norway
    24.0, # Boreal Baltic       — Estonia/Finland/Latvia/Lithuania/Sweden
    -5.0  # Iberian             — Portugal/Spain
  ),
  y = c(
    47.0,  # Central Alpine
    52.0,  # Atlantic
    51.5,  # Western Continental
    45.5,  # Eastern Alpine
    49.0,  # Eastern Continental
    38.0,  # Southeast Mediterranean
    47.2,  # Pannonian
    63.0,  # Scandinavian Highlands
    62.0,  # Boreal Baltic
    40.0   # Iberian
  ),
  stringsAsFactors = FALSE
)

bf_files    <- c(
  "Region_bf_equal.csv",
  "Region_bf_proportional.csv",
  "Region_bf_stratified.csv"
)
samp_labels <- c("equal", "proportional", "stratified")

for (i in seq_along(bf_files)) {
  
  bf_file <- bf_files[i]
  tag     <- samp_labels[i]
  message("Processing ", bf_file, " (", tag, ")")
  
  bf <- read_csv(bf_file, show_col_types = FALSE)
  
  arrows <- bf %>%
    filter(
      from %in% names(state_labels),
      to   %in% names(state_labels),
      from != to,
      bayes_factor > 100
    ) %>%
    group_by(from, to) %>%
    summarise(
      TransitionRate = mean(mean_rate, na.rm = TRUE),
      max_bf         = max(bayes_factor, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    mutate(rate_category = case_when(
      TransitionRate < 1                        ~ "<1",
      TransitionRate >= 1 & TransitionRate <= 2 ~ "1–2",
      TransitionRate > 2                        ~ ">3"
    )) %>%
    mutate(rate_category = factor(rate_category, levels = c("<1", "1–2", ">3")))
  
  if (nrow(arrows) == 0) {
    warning("No decisive transitions (BF > 100) for ", tag, " — skipping.")
    next
  }
  
  arrow_data <- arrows %>%
    left_join(cluster_centers, by = c("from" = "Region")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(cluster_centers, by = c("to" = "Region")) %>%
    rename(x_end = x, y_end = y)
  
  maps <- ggplot(europe_data) +
    geom_sf(aes(fill = Region), colour = "black", linewidth = 0.1) +
    scale_fill_manual(
      values   = region_pal,
      labels   = state_labels,
      name     = "HC Cluster",
      na.value = "grey90",
      drop     = FALSE
    ) +
    geom_curve(
      data = arrow_data,
      aes(
        x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        linewidth = rate_category
      ),
      curvature = 0.3,
      arrow     = arrow(angle = 15, type = "open", length = unit(0.12, "inches")),
      lineend   = "round",
      color     = "grey20"
    ) +
    scale_linewidth_manual(
      name   = "Transition Rate",
      values = c("<1" = 0.5, "1–2" = 1.0, ">3" = 2.0)
    ) +
    guides(fill = "none") +
    coord_sf(xlim = c(-12, 45), ylim = c(30, 72), expand = FALSE) +
    theme_classic() +
    theme(
      legend.position      = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background    = element_rect(fill = "white", colour = NA),
      axis.title           = element_blank(),
      axis.text            = element_blank(),
      axis.ticks           = element_blank(),
      panel.grid           = element_blank(),
      axis.line            = element_blank(),
      legend.title         = element_text(size = 20, face = "bold"),
      legend.text          = element_text(size = 18),
      legend.key.size      = unit(1.2, "cm")
    )
  
  out_png <- paste0("Region_HC_TransitionRates_", tag, ".png")
  ggsave(
    filename = out_png,
    plot     = maps,
    width    = 20,
    height   = 10,
    dpi      = 300,
    bg       = "transparent",
    units    = "in"
  )
  message("✓ ", out_png, " saved.")
}

