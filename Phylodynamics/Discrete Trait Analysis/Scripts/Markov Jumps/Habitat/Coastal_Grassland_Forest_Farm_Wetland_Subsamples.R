library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/Habitat")
options(scipen = 999)

jump_files <- c(
  "Habitat_Subsample1_jumpTimes.txt",
  "Habitat_Subsample2_jumpTimes.txt",
  "Habitat_Subsample3_jumpTimes.txt"
)

sub_labels <- c("Subsample1", "Subsample2", "Subsample3")

valid_states <- c("Coastal", "Wetland", "Farm", "Grassland")

src_cols <- list(
  "Coastal"   = c(scales::alpha("#0072B2", 0.25), "#0072B2"),
  "Farm"      = c(scales::alpha("#CD5C5C", 0.25), "#CD5C5C"),
  "Grassland" = c(scales::alpha("#26A69A", 0.25), "#26A69A"),
  "Wetland"   = c(scales::alpha("#006D2C", 0.25), "#006D2C")
)

strip_cols <- vapply(src_cols, function(x) x[2], character(1))

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

make_source_plot <- function(habitat_ave_complete, source_key, global_max,
                             show_xaxis = FALSE, show_labels = FALSE) {
  df <- habitat_ave_complete %>%
    filter(Source == source_key) %>%
    mutate(
      Sink_lab = factor(
        Sink_lab,
        levels = gsub("_", " ", valid_states),
        ordered = TRUE
      )
    )
  
  strip_col <- strip_cols[[as.character(source_key)]]
  
  ggplot(df, aes(x = factor(year), y = Sink_lab, fill = ave)) +
    geom_tile(color = "black", linewidth = 0.4) +
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
      axis.text.x       = if (show_xaxis) element_text(size = 16) else element_blank(),
      axis.ticks.x      = element_blank(),
      strip.text.y.left = element_text(face = "bold", size = 18, colour = strip_col)
    )
}

out_dir <- "Habitat"

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
  
  habitat_ave <- count_total %>%
    mutate(ave = n / state)
  
  habitat_ave_filtered <- habitat_ave %>%
    rename(Source = From, Sink = To) %>%
    filter(as.numeric(year) >= 2016, as.numeric(year) <= 2025) %>%
    mutate(year = as.integer(year))
  
  full_grid <- expand_grid(
    Source = valid_states,
    Sink   = valid_states,
    year   = 2016:2025
  ) %>% filter(Source != Sink)
  
  habitat_ave_complete <- full_grid %>%
    left_join(habitat_ave_filtered, by = c("Source", "Sink", "year")) %>%
    mutate(
      ave        = replace_na(ave, 0),
      Source     = factor(Source, levels = valid_states, ordered = TRUE),
      Sink       = factor(Sink,   levels = valid_states, ordered = TRUE),
      Source_lab = gsub("_", " ", Source),
      Sink_lab   = gsub("_", " ", Sink)
    )
  
  global_max <- max(habitat_ave_complete$ave, na.rm = TRUE)
  
  order_sources <- levels(habitat_ave_complete$Source)
  plots <- Map(
    function(h, idx) make_source_plot(
      habitat_ave_complete,
      source_key  = h,
      global_max  = global_max,
      show_xaxis  = idx == length(order_sources),
      show_labels = idx == length(order_sources)
    ),
    order_sources,
    seq_along(order_sources)
  )
  
  final_plot_habitat <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
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
  
  print(final_plot_habitat)
  
  png_file <- file.path(out_dir, paste0("CWFG_Transitions_tp_", tag, ".png"))
  pdf_file <- file.path(out_dir, paste0("CWFG_Transitions_tp_", tag, ".pdf"))
  
  ggsave(
    png_file,
    final_plot_habitat,
    bg     = "transparent",
    width  = 11,
    height = 8,
    dpi    = 600,
    units  = "in"
  )
  
  ggsave(
    pdf_file,
    final_plot_habitat,
    width       = 11,
    height      = 8,
    units       = "in",
    useDingbats = FALSE
  )
}
