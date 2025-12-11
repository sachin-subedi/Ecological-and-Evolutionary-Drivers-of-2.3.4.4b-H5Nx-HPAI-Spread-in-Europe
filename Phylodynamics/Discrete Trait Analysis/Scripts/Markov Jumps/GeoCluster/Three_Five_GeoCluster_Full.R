library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)
library(patchwork)
library(grid)

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

src_cols <- list(
  "GeoCluster_One"   = c(scales::alpha("#4DBBD5", 0.25), "#4DBBD5"),
  "GeoCluster_Two"   = c(scales::alpha("#F39B7F", 0.25), "#F39B7F"),
  "GeoCluster_Three" = c(scales::alpha("#3C5488", 0.25), "#3C5488"),
  "GeoCluster_Four"  = c(scales::alpha("#BCBD22", 0.25), "#BCBD22"),
  "GeoCluster_Five"  = c(scales::alpha("#1B9E77", 0.25), "#1B9E77")
)

strip_cols <- vapply(src_cols, function(x) x[2], character(1))

year_levels   <- as.character(2016:2025)
LEG_BARWIDTH  <- unit(80, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")

# Base theme
panel_theme <- theme_minimal(base_size = 18) +
  theme(
    panel.spacing.y    = unit(0, "lines"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(face = "bold", hjust = 1, size = 22),
    strip.background.y = element_blank(),
    axis.text.y        = element_text(size = 20),
    panel.grid         = element_blank(),
    plot.margin        = margin(0, 0, 0, 0)
  )

make_source_plot <- function(GeoCluster_ave_complete, source_key, global_max,
                             show_xaxis = FALSE, show_labels = FALSE) {
  
  df <- GeoCluster_ave_complete %>%
    filter(Source == source_key) %>%
    mutate(
      Sink_lab = factor(
        Sink_lab,
        levels = gsub("_", " ", valid_states),
        ordered = TRUE
      )
    )
  
  strip_col <- strip_cols[[source_key]]
  
  ggplot(df, aes(x = factor(year), y = Sink_lab, fill = ave)) +
    geom_tile(color = "black", linewidth = 0.6) +
    facet_grid(Source_lab ~ ., switch = "y", scales = "free_y", space = "free_y") +
    scale_x_discrete(drop = FALSE, limits = year_levels) +
    scale_fill_gradient(
      name   = if (source_key == "GeoCluster_Three") "Avg. Jumps per Year" else NULL,
      low    = src_cols[[source_key]][1],
      high   = src_cols[[source_key]][2],
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
      axis.text.x       = if (show_xaxis) element_text(size = 20) else element_blank(),
      axis.ticks.x      = element_blank(),
      strip.text.y.left = element_text(face = "bold", size = 22, colour = strip_col)
    )
}

out_dir <- "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/jumps/GeoCluster"

for (i in seq_along(jump_files)) {
  in_file <- jump_files[i]
  tag     <- sub_labels[i]
  
  message("Processing: ", in_file, " (", tag, ")")
  
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
  
  GeoCluster_ave <- count_total %>%
    mutate(ave = n / state)
  
  GeoCluster_ave_filtered <- GeoCluster_ave %>%
    rename(Source = From, Sink = To) %>%
    filter(as.numeric(year) >= 2016, as.numeric(year) <= 2025) %>%
    mutate(year = as.integer(year))
  
  full_grid <- expand_grid(
    Source = valid_states,
    Sink   = valid_states,
    year   = 2016:2025
  ) %>% filter(Source != Sink)
  
  GeoCluster_ave_complete <- full_grid %>%
    left_join(GeoCluster_ave_filtered, by = c("Source", "Sink", "year")) %>%
    mutate(
      ave        = replace_na(ave, 0),
      Source     = factor(Source, levels = valid_states, ordered = TRUE),
      Sink       = factor(Sink,   levels = valid_states, ordered = TRUE),
      Source_lab = gsub("_", " ", Source),
      Sink_lab   = gsub("_", " ", Sink)
    )
  
  global_max <- max(GeoCluster_ave_complete$ave, na.rm = TRUE)
  
  order_sources <- c("GeoCluster_Three", "GeoCluster_Five")
  
  plots <- Map(
    function(h, idx) make_source_plot(
      GeoCluster_ave_complete,
      source_key  = h,
      global_max  = global_max,
      show_xaxis  = idx == length(order_sources),
      show_labels = idx == length(order_sources)
    ),
    order_sources,
    seq_along(order_sources)
  )
  
  final_plot_GeoCluster <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
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
      legend.title         = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text          = element_text(size = 14),
      plot.margin          = margin(0, 0, 0, 0)
    )
  
  print(final_plot_GeoCluster)
  
  png_file <- file.path(out_dir, paste0("GeoCluster_ThreeFive_Transitions_tp_", tag, ".png"))
  pdf_file <- file.path(out_dir, paste0("GeoCluster_ThreeFive_Transitions_tp_", tag, ".pdf"))
  
  ggsave(
    png_file,
    final_plot_GeoCluster,
    bg     = "transparent",
    width  = 12,
    height = 12,
    dpi    = 600,
    units  = "in"
  )
  
  ggsave(
    pdf_file,
    final_plot_GeoCluster,
    width       = 12,
    height      = 12,
    units       = "in",
    useDingbats = FALSE
  )
}
