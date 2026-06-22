library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("results/jumps")
options(scipen = 999)

jump_files <- c(
  "Region_history_equal_combined.txt",
  "Region_history_proportional_combined.txt",
  "Region_history_stratified_combined.txt"
)

sub_labels <- c("equal", "proportional", "stratified")

source_states <- c(
  "HC_Cluster_1_Atlantic",
  "HC_Cluster_1_Continental",
  "HC_Cluster_2_Continental",
  "HC_Cluster_3_Boreal"
)

sink_states <- c(
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

valid_states <- union(source_states, sink_states)

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
  "HC_Cluster_1_Atlantic"    = c(scales::alpha("#CD5C5C", 0.25), "#CD5C5C"),
  "HC_Cluster_1_Continental" = c(scales::alpha("#BCBD22", 0.25), "#BCBD22"),
  "HC_Cluster_2_Continental"= c(scales::alpha("#0072B2",  0.25), "#0072B2"),
  "HC_Cluster_3_Boreal"      = c(scales::alpha("#006D2C", 0.25), "#006D2C")
)

strip_cols <- vapply(src_cols, function(x) x[2], character(1))
year_levels   <- as.character(2016:2025)
LEG_BARWIDTH  <- unit(80, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")

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

make_Region_plot <- function(Region_ave_complete, source_key, global_max,
                                 show_xaxis = FALSE, show_labels = FALSE) {
  df <- Region_ave_complete %>%
    filter(Source == source_key) %>%
    mutate(
      Sink_lab = factor(
        state_labels[as.character(Sink)],
        levels = state_labels[sink_states],
        ordered = TRUE
      )
    )
  
  strip_col <- strip_cols[[as.character(source_key)]]
  
  ggplot(df, aes(x = factor(year), y = Sink_lab, fill = ave)) +
    geom_tile(color = "black", linewidth = 0.2) +
    facet_grid(Source_lab ~ ., switch = "y", scales = "free_y", space = "free_y") +
    scale_x_discrete(drop = FALSE, limits = year_levels) +
    scale_fill_gradient(
      name   = if (source_key == valid_states[1]) "Avg. Jumps per Year" else NULL,
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
      axis.text.x  = if (show_xaxis) element_text(size = 20) else element_blank(),
      axis.ticks.x = element_blank(),
    )
}

out_dir <- "results/jumps/"

for (i in seq_along(jump_files)) {
  
  in_file <- jump_files[i]
  tag     <- sub_labels[i]
  
  message("Processing: ", in_file, " (", tag, ")")
  
  Region <- read.delim(in_file, sep = "", header = TRUE)
  
  Region <- Region %>%
    mutate(
      time    = 2025.2520547945205 - as.numeric(time),
      year    = format(date_decimal(time), "%Y"),
      from_to = paste(from, to, sep = ".")
    )
  
  state <- n_distinct(Region$state)
  
  count_total <- Region %>%
    group_by(year, from_to) %>%
    count(name = "n") %>%
    separate(from_to, into = c("From", "To"), sep = "\\.")
  
  count_total <- count_total %>%
    filter(From %in% source_states, To %in% sink_states)
  
  Region_ave <- count_total %>%
    mutate(ave = n / state)
  
  Region_ave_filtered <- Region_ave %>%
    rename(Source = From, Sink = To) %>%
    filter(as.numeric(year) >= 2016, as.numeric(year) <= 2025) %>%
    mutate(year = as.integer(year))
  
  full_grid <- expand_grid(
    Source = source_states,
    Sink   = sink_states,
    year   = 2016:2025
  ) %>% filter(Source != Sink)
  
  Region_ave_complete <- full_grid %>%
    left_join(Region_ave_filtered, by = c("Source", "Sink", "year")) %>%
    mutate(
      ave        = replace_na(ave, 0),
      Source     = factor(Source, levels = source_states, ordered = TRUE),
      Sink       = factor(Sink,   levels = sink_states,   ordered = TRUE),
      Source_lab = state_labels[as.character(Source)],
      Sink_lab   = state_labels[as.character(Sink)]
    )
  
  csv_file <- file.path(
    out_dir,
    paste0("Region_AvgJumpsPerYear_4states_", tag, ".csv")
  )
  
  Region_ave_complete %>%
    arrange(Source, Sink, year) %>%
    mutate(
      Source = as.character(Source),
      Sink   = as.character(Sink)
    ) %>%
    write.csv(csv_file, row.names = FALSE)
  
  global_max <- max(Region_ave_complete$ave, na.rm = TRUE)
  
  order_sources <- levels(Region_ave_complete$Source)
  plots <- Map(
    function(h, idx) make_Region_plot(
      Region_ave_complete,
      source_key  = h,
      global_max  = global_max,
      show_xaxis  = idx == length(order_sources),
      show_labels = idx == length(order_sources)
    ),
    order_sources,
    seq_along(order_sources)
  )
  
  final_plot_Region <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
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
  
  print(final_plot_Region)
  
  png_file <- file.path(out_dir, paste0("time_Region_heatmap_4states_", tag, ".png"))
  pdf_file <- file.path(out_dir, paste0("time_Region_heatmap_4states_", tag, ".pdf"))
  
  ggsave(
    filename = png_file,
    plot     = final_plot_Region,
    width    = 14, height = 16, dpi = 600, units = "in", type = "cairo"
  )
  
  ggsave(
    filename = pdf_file,
    plot     = final_plot_Region,
    width    = 14, height = 16, units = "in", useDingbats = FALSE
  )
}
