library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/jumps/GeoCluster")
options(scipen = 999)

jump_files <- c(
  "GeoCluster_Subsample1_jumpTimes.txt",
  "GeoCluster_Subsample2_jumpTimes.txt",
  "GeoCluster_Subsample3_jumpTimes.txt"
)
sub_labels <- c("Subsample1", "Subsample2", "Subsample3")

valid_states <- c(
  "GeoCluster_One",
  "GeoCluster_Two",
  "GeoCluster_Three",
  "GeoCluster_Four",
  "GeoCluster_Five"
)

keep_sources <- c("GeoCluster_Three", "GeoCluster_Five")

LABEL_MAP <- c(
  "GeoCluster_One"   = "South-Eastern",
  "GeoCluster_Two"   = "Central-Eastern",
  "GeoCluster_Three" = "Central-Western",
  "GeoCluster_Four"  = "Northern",
  "GeoCluster_Five"  = "Central-Southern"
)

src_cols <- list(
  "GeoCluster_One"   = c(scales::alpha("#4DBBD5", 0.25), "#4DBBD5"),
  "GeoCluster_Two"   = c(scales::alpha("#F39B7F", 0.25), "#F39B7F"),
  "GeoCluster_Three" = c(scales::alpha("#3C5488", 0.25), "#3C5488"),
  "GeoCluster_Four"  = c(scales::alpha("#BCBD22", 0.25), "#BCBD22"),
  "GeoCluster_Five"  = c(scales::alpha("#1B9E77", 0.25), "#1B9E77")
)

strip_cols <- c(
  "GeoCluster_One"   = "#4DBBD5",
  "GeoCluster_Two"   = "#F39B7F",
  "GeoCluster_Three" = "#3C5488",
  "GeoCluster_Four"  = "#BCBD22",
  "GeoCluster_Five"  = "#1B9E77"
)

year_levels   <- as.character(2016:2025)
LEG_BARWIDTH  <- unit(80, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")

panel_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.spacing.y    = unit(0, "lines"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(face = "bold", hjust = 1, size = 22),
    strip.background.y = element_blank(),
    axis.text.y        = element_text(size = 16),
    panel.grid         = element_blank(),
    plot.margin        = margin(0, 0, 0, 0)
  )

make_source_plot <- function(geo_ave_complete, source_key, global_max,
                             show_xaxis = FALSE, show_labels = FALSE) {
  df <- geo_ave_complete %>%
    filter(Source == source_key) %>%
    mutate(
      Sink_lab = factor(
        Sink_lab,
        levels  = unname(LABEL_MAP[valid_states]),
        ordered = TRUE
      )
    )
  
  strip_col <- strip_cols[[as.character(source_key)]]
  
  ggplot(df, aes(x = factor(year), y = Sink_lab, fill = ave)) +
    geom_tile(color = "black", linewidth = 0.4) +
    facet_grid(Source_lab ~ ., switch = "y", scales = "free_y", space = "free_y") +
    scale_x_discrete(drop = FALSE, limits = year_levels) +
    scale_fill_gradient(
      name   = if (source_key == keep_sources[1]) "Avg. Jumps per Year" else NULL,
      low    = src_cols[[as.character(source_key)]][1],
      high   = src_cols[[as.character(source_key)]][2],
      limits = c(0, global_max),
      guide  = guide_colorbar(
        direction      = "horizontal",
        title.position = "top",
        barwidth       = LEG_BARWIDTH,
        barheight      = LEG_BARHEIGHT,
        ticks          = show_labels,
        label          = show_labels
      )
    ) +
    labs(x = if (show_xaxis) "" else NULL, y = NULL) +
    panel_theme +
    theme(
      axis.text.x       = if (show_xaxis) element_text(size = 16) else element_blank(),
      axis.ticks.x      = element_blank(),
      strip.text.y.left = element_text(face = "bold", size = 18, colour = strip_col)
    )
}

out_dir <- "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/jumps/GeoCluster/Regions"

for (i in seq_along(jump_files)) {
  in_file <- jump_files[i]
  tag     <- sub_labels[i]
  
  jumps <- read.delim(in_file, sep = "", header = TRUE)
  
  jumps <- jumps %>%
    mutate(
      time    = 2025.2520547945205 - as.numeric(time),
      year    = format(date_decimal(time), "%Y"),
      from_to = paste(from, to, sep = ".")
    )
  
  state <- n_distinct(jumps$state)
  
  count_total <- jumps %>%
    group_by(year, from_to) %>%
    count(name = "n") %>%
    separate(from_to, into = c("From", "To"), sep = "\\.")
  
  count_total <- count_total %>%
    filter(From %in% valid_states, To %in% valid_states)
  
  geo_ave <- count_total %>%
    mutate(ave = n / state)
  
  geo_ave_filtered <- geo_ave %>%
    rename(Source = From, Sink = To) %>%
    filter(as.numeric(year) >= 2016, as.numeric(year) <= 2025) %>%
    mutate(year = as.integer(year)) %>%
    filter(Source %in% keep_sources)
  
  full_grid <- expand_grid(
    Source = keep_sources,
    Sink   = valid_states,
    year   = 2016:2025
  ) %>% filter(Source != Sink)
  
  geo_ave_complete <- full_grid %>%
    left_join(geo_ave_filtered, by = c("Source", "Sink", "year")) %>%
    mutate(
      ave    = replace_na(ave, 0),
      Source = factor(Source, levels = keep_sources, ordered = TRUE),
      Sink   = factor(Sink,   levels = valid_states, ordered = TRUE),
      Source_lab = factor(
        unname(LABEL_MAP[as.character(Source)]),
        levels  = unname(LABEL_MAP[keep_sources]),
        ordered = TRUE
      ),
      Sink_lab = factor(
        unname(LABEL_MAP[as.character(Sink)]),
        levels  = unname(LABEL_MAP[valid_states]),
        ordered = TRUE
      )
    )
  
  csv_file <- file.path(out_dir, paste0("GeoCluster_AvgJumpsPerYear_", tag, "_Src3_5.csv"))
  
  geo_ave_complete %>%
    arrange(Source, Sink, year) %>%
    mutate(
      Source = as.character(Source),
      Sink   = as.character(Sink)
    ) %>%
    write.csv(csv_file, row.names = FALSE)
  
  global_max <- max(geo_ave_complete$ave, na.rm = TRUE)
  order_sources <- levels(geo_ave_complete$Source)
  
  plots <- Map(
    function(h, idx) make_source_plot(
      geo_ave_complete,
      source_key  = h,
      global_max  = global_max,
      show_xaxis  = idx == length(order_sources),
      show_labels = idx == length(order_sources)
    ),
    order_sources,
    seq_along(order_sources)
  )
  
  final_plot_geo <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
    plot_layout(guides = "collect") &
    theme(
      legend.position      = "bottom",
      legend.box           = "vertical",
      legend.direction     = "vertical",
      legend.box.just      = "center",
      legend.justification = "center",
      legend.spacing.y     = unit(1, "pt"),
      legend.margin        = margin(0, 0, 0, 0),
      legend.background    = element_rect(fill = "white", colour = NA),
      legend.title         = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.text          = element_text(size = 12),
      plot.margin          = margin(0, 0, 0, 0)
    )
  
  print(final_plot_geo)
  
  png_file <- file.path(out_dir, paste0("Full_Region_Transitions_tp_", tag, "_Src3_5.png"))
  pdf_file <- file.path(out_dir, paste0("Full_Region_Transitions_tp_", tag, "_Src3_5.pdf"))
  
  ggsave(
    png_file,
    final_plot_geo,
    bg     = "transparent",
    width  = 8,
    height = 14,
    dpi    = 600,
    units  = "in"
  )
  
  ggsave(
    pdf_file,
    final_plot_geo,
    width       = 8,
    height      = 14,
    units       = "in",
    useDingbats = FALSE
  )
}
