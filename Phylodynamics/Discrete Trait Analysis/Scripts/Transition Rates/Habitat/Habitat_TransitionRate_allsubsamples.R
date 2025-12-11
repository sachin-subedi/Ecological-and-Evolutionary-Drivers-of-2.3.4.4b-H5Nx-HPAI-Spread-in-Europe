options(scipen = 999)

library(dplyr)
library(ggplot2)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/rates/Habitat")

bf_files <- c(
  "Habitat_bf_Subsample1.csv",
  "Habitat_bf_Subsample2.csv",
  "Habitat_bf_Subsample3.csv"
)
tags <- c("Subsample1", "Subsample2", "Subsample3")

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
      From = from,
      To   = to,
      Rate = mean_rate
    )
  
  df_bar <- df %>%
    dplyr::semi_join(supported_pairs, by = c("From", "To")) %>%
    dplyr::filter(bayes_factor >= 10) %>% 
    dplyr::group_by(From, To) %>%
    dplyr::summarise(Rate = mean(Rate), .groups = "drop") %>%
    dplyr::mutate(
      Transition = paste(From, "\u2192", To)
    ) %>%
    dplyr::arrange(desc(Rate))
  
  df_bar$Transition <- factor(df_bar$Transition,
                              levels = rev(df_bar$Transition))
  
  p_bar <- ggplot(df_bar,
                  aes(x = Transition, y = Rate, fill = From)) +
    geom_col() +
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
  
  out_name <- paste0("Habitat_Rates_barplot_", tag_i, "_light0.5.png")
  ggsave(
    filename = out_name,
    plot     = p_bar,
    width    = 8,
    height   = 10,
    bg       = "white",
    dpi      = 300
  )
  
  message("✓ saved: ", out_name)
}

