library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/jumps/Host")
options(scipen = 999)

# ── Subsample files ----------------------------------------------------------
jump_files <- c(
  "Host_Subsample1_jumpTimes.txt",
  "Host_Subsample2_jumpTimes.txt",
  "Host_Subsample3_jumpTimes.txt"
)

sub_labels <- c("Subsample1", "Subsample2", "Subsample3")

# ── Host states & colors -----------------------------------------------------
valid_states <- c("Domestic_Bird", "Domestic_Mammal", "Human", "Wild_Bird", "Wild_Mammal")

src_cols <- list(
  "Domestic_Bird"   = c(scales::alpha("#CD5C5C", 0.25), "#CD5C5C"),
  "Domestic_Mammal" = c(scales::alpha("#6A3D9A", 0.25), "#6A3D9A"),
  "Human"           = c(scales::alpha("#1B9E77", 0.25), "#1B9E77"),
  "Wild_Bird"       = c(scales::alpha("#3C5488", 0.25), "#3C5488"),
  "Wild_Mammal"     = c(scales::alpha("#BCBD22", 0.25), "#BCBD22")
)

# Use the "high" color of each palette for strip text coloring
strip_cols <- vapply(src_cols, function(x) x[2], character(1))

year_levels   <- as.character(2016:2025)
LEG_BARWIDTH  <- unit(80, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")

# Base theme
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

# Helper: build one plot per Source host
make_source_plot <- function(host_ave_complete, source_key, global_max,
                             show_xaxis = FALSE, show_labels = FALSE) {
  df <- host_ave_complete %>%
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
      # strip.background.y = element_rect(fill = scales::alpha(strip_col, 0.10), colour = NA)
    )
}

# ── Output directory ---------------------------------------------------------
out_dir <- "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/jumps/Host"

# ── Loop over the 3 subsamples ----------------------------------------------
for (i in seq_along(jump_files)) {
  in_file <- jump_files[i]
  tag     <- sub_labels[i]
  
  message("Processing: ", in_file, " (", tag, ")")
  
  # Step 0: Read data
  jumps <- read.delim(in_file, sep = "", header = TRUE)
  
  # Step 1: Add time & year
  jumps <- jumps %>%
    mutate(
      time    = 2025.2520547945205 - as.numeric(time),
      year    = format(date_decimal(time), "%Y"),
      from_to = paste(from, to, sep = ".")
    )
  
  # Step 2: number of unique states (used for averaging)
  state <- n_distinct(jumps$state)
  
  # Step 3: Summarize transitions per year
  count_total <- jumps %>%
    group_by(year, from_to) %>%
    count(name = "n") %>%
    separate(from_to, into = c("From", "To"), sep = "\\.")
  
  # Step 4: keep only Host states of interest
  count_total <- count_total %>%
    filter(From %in% valid_states, To %in% valid_states)
  
  # Step 5: Average transitions
  host_ave <- count_total %>%
    mutate(ave = n / state)
  
  host_ave_filtered <- host_ave %>%
    rename(Source = From, Sink = To) %>%
    filter(as.numeric(year) >= 2016, as.numeric(year) <= 2025) %>%
    mutate(year = as.integer(year))
  
  # Step 6: Full grid (no self transitions)
  full_grid <- expand_grid(
    Source = valid_states,
    Sink   = valid_states,
    year   = 2016:2025
  ) %>% filter(Source != Sink)
  
  # Step 7: Complete & label-tidy
  host_ave_complete <- full_grid %>%
    left_join(host_ave_filtered, by = c("Source", "Sink", "year")) %>%
    mutate(
      ave        = replace_na(ave, 0),
      Source     = factor(Source, levels = valid_states, ordered = TRUE),
      Sink       = factor(Sink,   levels = valid_states, ordered = TRUE),
      Source_lab = gsub("_", " ", Source),
      Sink_lab   = gsub("_", " ", Sink)
    )
  
  # Per-subsample max
  global_max <- max(host_ave_complete$ave, na.rm = TRUE)
  
  # Build all panels (only bottom shows x-axis & legend labels)
  order_sources <- levels(host_ave_complete$Source)
  plots <- Map(
    function(h, idx) make_source_plot(
      host_ave_complete,
      source_key  = h,
      global_max  = global_max,
      show_xaxis  = idx == length(order_sources),
      show_labels = idx == length(order_sources)
    ),
    order_sources,
    seq_along(order_sources)
  )
  
  final_plot_host <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
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
  
  print(final_plot_host)
  
  # Save with subsample tag
  png_file <- file.path(out_dir, paste0("Host_Transitions_tp_", tag, ".png"))
  pdf_file <- file.path(out_dir, paste0("Host_Transitions_tp_", tag, ".pdf"))
  
  ggsave(
    png_file,
    final_plot_host,
    bg     = "transparent",
    width  = 8,
    height = 11,
    dpi    = 600,
    units  = "in"
  )
  
  ggsave(
    pdf_file,
    final_plot_host,
    width       = 8,
    height      = 11,
    units       = "in",
    useDingbats = FALSE
  )
}
