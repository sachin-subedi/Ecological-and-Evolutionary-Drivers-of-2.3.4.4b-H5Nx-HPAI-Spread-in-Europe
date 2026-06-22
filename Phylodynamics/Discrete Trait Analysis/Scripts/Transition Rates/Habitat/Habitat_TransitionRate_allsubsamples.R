options(scipen = 999)
library(dplyr)
library(ggplot2)

setwd("Habitat/")

bf_files <- c(
  "Habitat_bf_CI_equal.csv",
  "Habitat_bf_CI_proportional.csv",
  "Habitat_bf_CI_stratified.csv"
)
tags <- c("equal", "proportional", "stratified")

habitat_cols_solid <- c(
  "Coastal"        = "#0072B2",
  "Farm"           = "#CD5C5C",
  "Forest"         = "#BCBD22",
  "Grassland"      = "#26A69A",
  "Human_Modified" = "#56B4E9",
  "Marine"         = "#CC79A7",
  "Shrubland"      = "lightslategray",
  "Urban"          = "#E6842A",
  "Woodland"       = "#6B4C1B",
  "Rock"           = "#7570B3",
  "Wetland"        = "#006D2C"
)

lighten_palette <- function(cols, factor = 0.5) {
  sapply(cols, function(color) {
    rgb_col <- col2rgb(color)
    new_col <- rgb_col + (255 - rgb_col) * factor
    rgb(new_col[1], new_col[2], new_col[3], maxColorValue = 255)
  })
}
habitat_cols_light <- lighten_palette(habitat_cols_solid, factor = 0.5)

supported_pairs <- tibble::tribble(
  ~From,       ~To,
  "Wetland",   "Forest",
  "Wetland",   "Coastal",
  "Wetland",   "Farm",
  "Wetland",   "Grassland",
  "Grassland", "Coastal",
  "Grassland", "Farm",
  "Farm",      "Wetland",
  "Coastal",   "Forest",
  "Coastal",   "Grassland",
  "Coastal",   "Marine"
)

for (i in seq_along(bf_files)) {
  
  file_i <- bf_files[i]
  tag_i  <- tags[i]
  
  message("Processing: ", file_i)
  
  df <- read.csv(file_i) %>%
    rename(
      From  = from,
      To    = to,
      Rate  = mean_real_rate, 
      Lower = ci_lower,
      Upper = ci_upper
    )
  
  df_bar <- df %>%
    dplyr::semi_join(supported_pairs, by = c("From", "To")) %>%
    dplyr::filter(bayes_factor >= 10) %>%
    dplyr::group_by(From, To) %>%
    dplyr::summarise(
      Rate  = mean(Rate),
      Lower = mean(Lower),
      Upper = mean(Upper),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Transition = paste(From, "\u2192", To)
    ) %>%
    dplyr::arrange(desc(Rate))
  
  df_bar$Transition <- factor(df_bar$Transition,
                              levels = rev(df_bar$Transition))
  
  p_bar <- ggplot(df_bar,
                  aes(x = Transition, y = Rate, fill = From)) +
    geom_col() +
    geom_errorbar(                              # ← 95% credible interval error bars
      aes(ymin = Lower, ymax = Upper),
      width     = 0.3,
      linewidth = 0.7,
      color     = "grey30"
    ) +
    coord_flip() +
    scale_fill_manual(values = habitat_cols_light, name = "Source habitat") +
    labs(
      x = "Habitat Transitions",
      y = "Transition rate"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 16),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold")
    )
  
  print(p_bar)
  
  out_name <- paste0("Habitat_Rates_barplot_", tag_i, "_light0.5_CI.png")
  ggsave(
    filename = out_name,
    plot     = p_bar,
    width    = 8,
    height   = 10,
    bg       = "white",
    dpi      = 300
  ) 
  
  message("\u2713 saved: ", out_name)
}

df %>%
  mutate(
    hpd_width = hpd_upper - hpd_lower,
    ci_width  = Upper - Lower,
    diff      = ci_width - hpd_width
  ) %>%
  select(From, To, hpd_lower, hpd_upper, Lower, Upper, hpd_width, ci_width, diff) %>%
  arrange(desc(abs(diff)))
