library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(RColorBrewer)
library(grid)
library(lwgeom)
library(cowplot)

setwd("/GeoCluster")
options(scipen = 999)

lighten_color <- function(color, factor = 0.6) {
  col <- col2rgb(color)
  col <- col + (255 - col) * factor
  rgb(t(col), maxColorValue = 255)
}

zone <- c("GeoCluster_One", "GeoCluster_Two", "GeoCluster_Three",
          "GeoCluster_Four", "GeoCluster_Five")

countries <- list(
  c("Serbia", "Romania", "North Macedonia", "Kosovo", "Greece", "Cyprus",
    "Bulgaria", "Albania", "Montenegro"),                                            # GeoCluster_One
  c("Ukraine", "Slovakia", "Poland", "Moldova", "Hungary"),                          # GeoCluster_Two
  c("United Kingdom", "Switzerland", "Spain", "Portugal", "Netherlands", "Luxembourg",
    "Ireland", "Iceland", "Germany", "France", "Denmark", "Belgium"),                # GeoCluster_Three
  c("Sweden", "Norway", "Lithuania", "Latvia", "Finland", "Estonia", "Belarus"),    # GeoCluster_Four
  c("Slovenia", "Italy", "Czechia", "Croatia", "Bosnia and Herz.", "Austria")       # GeoCluster_Five
)

country_data <- data.frame(
  Country = unlist(countries),
  zone    = rep(zone, sapply(countries, length))
)

europe <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(name != "Russia")

europe_data <- europe %>%
  left_join(country_data, by = c("name" = "Country"))

europe_data$zone <- factor(
  europe_data$zone,
  levels = zone,
  labels = zone
)

region_pal <- c(
  "GeoCluster_One"   = "#4DBBD5",
  "GeoCluster_Two"   = "#F39B7F",
  "GeoCluster_Three" = "#3C5488",
  "GeoCluster_Four"  = "#BCBD22",
  "GeoCluster_Five"  = "#1B9E77"
)

region_pal_light <- sapply(region_pal, lighten_color, factor = 0.2)

zone_centers <- data.frame(
  zone = zone,
  x = c(25.7, 23.0, -1.5, 10.8, 12.5),
  y = c(45.9, 49.5, 52.5, 60.5, 45.5)
)

bf_files <- c(
  "GeoCluster_bf_Subsample1.csv",
  "GeoCluster_bf_Subsample2.csv",
  "GeoCluster_bf_Subsample3.csv"
)

sub_labels <- c("Subsample1", "Subsample2", "Subsample3")

for (i in seq_along(bf_files)) {
  
  bf_file <- bf_files[i]
  tag     <- sub_labels[i]
  message("Processing ", bf_file, " (", tag, ")")
  
  bf <- read.csv(bf_file, header = TRUE)
  
  arrows <- bf %>%
    filter(
      from %in% zone,
      to   %in% zone,
      from != to,
      bayes_factor > 100
    ) %>%
    group_by(from, to) %>%
    summarise(
      TransitionRate = mean(mean_rate, na.rm = TRUE),
      max_bf         = max(bayes_factor, na.rm = TRUE),
      .groups        = "drop"
    )

  if (nrow(arrows) == 0) {
    warning("No decisive transitions (BF > 100) for ", tag, " — skipping.")
    next
  }

  arrows <- arrows %>%
    mutate(rate_category = dplyr::case_when(
      TransitionRate < 1                        ~ "<1",
      TransitionRate >= 1 & TransitionRate <= 2 ~ "1–2",
      TransitionRate > 2                        ~ ">3"
    )) %>%
    mutate(rate_category = factor(rate_category, levels = c("<1", "1–2", ">3")))
  
  arrow_data <- arrows %>%
    left_join(zone_centers, by = c("from" = "zone")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(zone_centers, by = c("to" = "zone")) %>%
    rename(x_end = x, y_end = y) %>%
    tidyr::fill(c(x_start, y_start, x_end, y_end), .direction = "down")
  
  maps <- ggplot(europe_data) +
    geom_sf(aes(fill = zone), colour = "black", size = 0.1) +
    scale_fill_manual(values = region_pal_light, guide = "none") +
    
    geom_curve(
      data = arrow_data,
      aes(x = x_start, y = y_start,
          xend = x_end, yend = y_end,
          size = rate_category),
      curvature = 0.5,
      arrow = arrow(angle = 15, type = "open", length = unit(0.12, "inches")),
      lineend = "round",
      color = "grey20"
    ) +
    scale_size_manual(
      name = "Transition Rate",
      values = c("<1" = 0.5, "1–2" = 1.0, ">3" = 2)
    ) +
    coord_sf(xlim = c(-10, 60), ylim = c(30, 82), expand = FALSE) +
    
    theme_classic() +
    theme(
      legend.position       = c(0.0, 0.75),
      legend.justification  = c(0, 0.5),
      legend.background     = element_rect(fill = "white", colour = NA),
      axis.title            = element_blank(),
      axis.text             = element_blank(),
      axis.ticks            = element_blank(),
      panel.grid            = element_blank(),
      axis.line             = element_blank(),
      legend.title          = element_text(size = 18),
      legend.text           = element_text(size = 16)
    )
  
  print(maps)
  
  out_png <- paste0("GeoCluster_TransitionRates_", tag, ".png")
  ggsave(
    filename = out_png,
    plot = maps,
    width = 16,
    height = 10,
    dpi = 300,
    bg = "transparent",
    units = "in"
  )
}
