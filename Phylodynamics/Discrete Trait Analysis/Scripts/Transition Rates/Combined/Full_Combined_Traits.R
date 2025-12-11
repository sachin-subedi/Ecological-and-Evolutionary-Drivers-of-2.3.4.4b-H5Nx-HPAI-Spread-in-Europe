options(scipen = 999)

library(dplyr)
library(ggplot2)
library(ggnewscale)
library(scales)
library(grid)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/rates/Combined/")

input_files <- c(
  "HG_bf_Subsample1.csv",
  "HG_bf_Subsample2.csv",
  "HG_bf_Subsample3.csv"
)
tags <- c("Subsample1", "Subsample2", "Subsample3")

rename_labels <- function(x) {
  x_chr <- as.character(x)
  x_chr[grepl("^GAP", x_chr) | is.na(x_chr)] <- ""
  x_chr
}

hab_codes <- c("Coastal", "Wetland", "Farm", "Grassland", "Forest")
suffixes  <- 1:5

nodes <- unlist(
  lapply(suffixes, function(s) {
    group <- paste0(hab_codes, s)
    if (s < max(suffixes)) c(group, paste0("GAP", s)) else group
  })
)

suffix_cols <- list(
  "1" = c(alpha("#4DBBD5", 0.25), "#4DBBD5"),
  "2" = c(alpha("#F39B7F", 0.25), "#F39B7F"),
  "3" = c(alpha("#3C5488", 0.25), "#3C5488"),
  "4" = c(alpha("#BCBD22", 0.25), "#BCBD22"),
  "5" = c(alpha("#1B9E77", 0.25), "#1B9E77")
)

global_max   <- 3
LEG_BARWIDTH  <- unit(80, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")

panel_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.grid  = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.position      = "bottom",
    legend.box           = "vertical",
    legend.direction     = "horizontal",
    legend.box.just      = "center",
    legend.justification = "center",
    legend.spacing.y     = unit(1, "pt"),
    legend.margin        = margin(0, 0, 0, 0),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.title         = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text          = element_text(size = 16)
  )

for (k in seq_along(input_files)) {
  
  file_in <- input_files[k]
  tag     <- tags[k]
  message("Processing: ", file_in, " (", tag, ")")
  
  df <- read.csv(file_in)
  
  df <- df %>%
    rename(
      From = from,
      To   = to,
      Rate = mean_rate
    ) %>%
    mutate(
      Rate_plot = dplyr::case_when(
        From == To ~ NA_real_,
        Rate < 1   ~ 0,
        Rate < 2   ~ 1,
        Rate < 3   ~ 2,
        TRUE       ~ 3
      ),
      
      BF_cat = dplyr::case_when(
        From != To &
          !grepl("^GAP", From) &
          !grepl("^GAP", To) &
          bayes_factor >= 10 & bayes_factor < 30  ~ "10–30 (strong)",
        From != To &
          !grepl("^GAP", From) &
          !grepl("^GAP", To) &
          bayes_factor >= 30 & bayes_factor < 100 ~ "30–100 (very strong)",
        From != To &
          !grepl("^GAP", From) &
          !grepl("^GAP", To) &
          bayes_factor >= 100                     ~ ">100 (decisive)",
        TRUE ~ NA_character_
      ),
      BF_star = dplyr::case_when(
        BF_cat == "10–30 (strong)"          ~ "*",
        BF_cat == "30–100 (very strong)"    ~ "**",
        BF_cat == ">100 (decisive)"         ~ "***",
        TRUE ~ ""
      )
    )
  
  df$From <- factor(df$From, levels = nodes)
  df$To   <- factor(df$To,   levels = nodes)
  
  df_within <- df %>%
    filter(
      !grepl("^GAP", From),
      !grepl("^GAP", To)
    ) %>%
    mutate(
      suf_from = sub(".*?(\\d)$", "\\1", as.character(From)),
      suf_to   = sub(".*?(\\d)$", "\\1", as.character(To))
    ) %>%
    filter(suf_from == suf_to)

  p <- ggplot() +
    scale_x_discrete(
      drop   = FALSE,
      labels = rename_labels,
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      drop   = FALSE,
      labels = rename_labels,
      expand = c(0, 0)
    ) +
    labs(x = "Sink", y = "Source") +
    panel_theme +
    theme(
      axis.text.x  = element_text(size = 20, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 20),
      axis.title.x = element_text(size = 22, face = "bold"),  # Sink
      axis.title.y = element_text(size = 22, face = "bold")   # Source
    )
  
  for (s in as.character(1:5)) {
    subdat <- df %>%
      filter(sub(".*?(\\d)$", "\\1", as.character(From)) == s)
    
    ord <- as.integer(s)
    
    p <- p +
      geom_tile(
        data = subdat,
        aes(x = To, y = From, fill = Rate_plot),
        color = "black",
        linewidth = 0.5
      ) +
      scale_fill_gradient(
        name   = if (s == "1") "Transition rate" else NULL,  # title only once
        low    = suffix_cols[[s]][1],
        high   = suffix_cols[[s]][2],
        limits = c(0, global_max),
        breaks = c(0, 1, 2, 3),   # show 0,1,2,3 in legend
        na.value = "white",
        guide  = guide_colorbar(
          order         = ord,
          direction      = "horizontal",
          title.position = "top",
          barwidth       = LEG_BARWIDTH,
          barheight      = LEG_BARHEIGHT,
          ticks          = (s == "5"),
          label          = (s == "5"),
          label.position = "bottom"
        )
      )
    
    if (s != "5") p <- p + ggnewscale::new_scale_fill()
  }
  p <- p +
    geom_tile(
      data = df_within,
      aes(x = To, y = From),
      fill = NA,
      colour = "black",
      linewidth = 0.4
    )
  
  p <- p +
    geom_text(
      data = df %>% filter(BF_star != ""),
      aes(x = To, y = From, label = BF_star),
      size = 6,
      fontface = "bold"
    )
  
  bf_legend_df <- tibble(
    To     = factor(levels(df$To)[1],   levels = levels(df$To)),
    From   = factor(levels(df$From)[1], levels = levels(df$From)),
    BF_cat = factor(
      c("10–30 (strong)",
        "30–100 (very strong)",
        ">100 (decisive)"),
      levels = c("10–30 (strong)",
                 "30–100 (very strong)",
                 ">100 (decisive)")
    )
  )
  
  p <- p +
    geom_text(
      data  = bf_legend_df,
      aes(x = To, y = From, colour = BF_cat, label = "a"),
      size  = 6,
      fontface = "bold",
      alpha = 0,
      show.legend = TRUE
    ) +
    scale_colour_manual(
      name   = "Bayes factor",
      values = c("black", "black", "black"),
      guide  = guide_legend(
        order          = 6,
        title.position = "top",
        direction      = "vertical",
        nrow           = 3,
        byrow          = TRUE,
        override.aes   = list(
          alpha = 1,
          label = c("*", "**", "***")
        )
      )
    )
  
  print(p)
  
  out_file <- paste0(
    "HG_Rates_Full_heatmap_allSubsamples_",
    tag, ".png"
  )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 18,
    height = 14,
    bg = "white",
    dpi = 300
  )
}

