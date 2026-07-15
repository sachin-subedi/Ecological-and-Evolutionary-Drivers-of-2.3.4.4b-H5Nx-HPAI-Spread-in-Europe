options(scipen = 999)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(scales)
library(grid)
datasets <- list(
  list(input_file = "HG_bf_equal.csv",        tag = "equal"),
  list(input_file = "HG_bf_proportional.csv", tag = "proportional"),
  list(input_file = "HG_bf_stratified.csv",   tag = "stratified")
)
hc_region_states <- list(
  "HC1_Alp" = c("CM", "FO", "GW", "WT"),
  "HC1_Atl" = c("CM", "FA", "FO", "GW", "UB", "WT"),
  "HC1_Con" = c("CM", "FA", "FO", "GW", "WT"),
  "HC2_Alp" = c("CM", "FA", "FO", "GW", "WT"),
  "HC2_Con" = c("CM", "FA", "FO", "GW", "UB", "WT"),
  "HC2_Med" = c("CM", "FA", "FO", "GW", "UB", "WT"),
  "HC2_Pan" = c("CM", "FA", "FO", "GW", "WT"),
  "HC3_Alp" = c("CM", "FA", "FO", "GW", "WT"),
  "HC3_Bor" = c("CM", "FA", "FO", "GW", "WT"),
  "HC4_Med" = c("CM", "FA", "FO", "GW", "WT")
)

group_cols <- list(
  "HC1_Alp" = c(alpha("#0072B2",        0.25), "#0072B2"),
  "HC1_Atl" = c(alpha("#CD5C5C",        0.25), "#CD5C5C"),
  "HC1_Con" = c(alpha("#BCBD22",        0.25), "#BCBD22"),
  "HC2_Alp" = c(alpha("#26A69A",        0.25), "#26A69A"),
  "HC2_Con" = c(alpha("#CC79A7",        0.25), "#CC79A7"),
  "HC2_Med" = c(alpha("lightslategray", 0.25), "lightslategray"),
  "HC2_Pan" = c(alpha("#6B4C1B",        0.25), "#6B4C1B"),
  "HC3_Alp" = c(alpha("#7570B3",        0.25), "#7570B3"),
  "HC3_Bor" = c(alpha("#006D2C",        0.25), "#006D2C"),
  "HC4_Med" = c(alpha("#56B4E9",        0.25), "#56B4E9")
)

group_titles <- c(
  "HC1_Alp" = "Central\nAlpine",
  "HC1_Atl" = "Atlantic",
  "HC1_Con" = "Western\nContinental",
  "HC2_Alp" = "Eastern\nAlpine",
  "HC2_Con" = "Eastern\nContinental",
  "HC2_Med" = "Southeast\nMediterranean",
  "HC2_Pan" = "Pannonian",
  "HC3_Alp" = "Scandinavian\nHighlands",
  "HC3_Bor" = "Boreal\nBaltic",
  "HC4_Med" = "Iberian"
)
groups_keep <- names(hc_region_states)
rename_labels <- function(x) {
  hab_map <- c(CM = "Coastal", FA = "Farm", FO = "Forest",
               GW = "Grassland", WT = "Wetland", UB = "Urban")
  x_chr <- as.character(x)
  x_chr[grepl("^GAP", x_chr) | is.na(x_chr)] <- ""
  idx <- x_chr != "" & !grepl("^GAP", x_chr)
  if (any(idx)) {
    codes        <- sub("^HC\\d+_[A-Za-z]+_", "", x_chr[idx])
    x_chr[idx]  <- ifelse(codes %in% names(hab_map), hab_map[codes], codes)
  }
  x_chr
}
get_group <- function(x) sub("^(HC\\d+_[A-Za-z]+)_.*$", "\\1", as.character(x))
global_max    <- 3
LEG_BARWIDTH  <- unit(55, "pt")
LEG_BARHEIGHT <- unit(6,  "pt")
panel_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.grid  = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 5),
    legend.position      = "bottom",
    legend.box           = "vertical",
    legend.direction     = "horizontal",
    legend.box.just      = "center",
    legend.justification = "center",
    legend.spacing.y     = unit(1, "pt"),
    legend.margin        = margin(0, 0, 0, 0),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.title         = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text          = element_text(size = 14)
  )

nodes <- unlist(
  lapply(seq_along(groups_keep), function(i) {
    g   <- groups_keep[i]
    grp <- paste0(g, "_", hc_region_states[[g]])
    if (i < length(groups_keep)) c(grp, paste0("GAP", i)) else grp
  })
)
group_y_mids <- sapply(groups_keep, function(g) {
  pos <- match(paste0(g, "_", hc_region_states[[g]]), nodes)
  mean(c(min(pos), max(pos)))
})
side_label_df <- tibble(
  label = unname(group_titles[groups_keep]),
  ymid  = unname(group_y_mids)
)
for (ds in datasets) {
  
  input_file <- ds$input_file
  tag        <- ds$tag
  
  message("Processing: ", input_file, " (", tag, ")")
  
  df <- read.csv(input_file)
  
  df <- df %>%
    rename(From = from, To = to, Rate = median_rate) %>%
    mutate(
      Rate_plot = dplyr::case_when(
        From == To ~ NA_real_,
        Rate < 1   ~ 0,
        Rate < 2   ~ 1,
        Rate < 3   ~ 2,
        TRUE       ~ 3
      ),
      BF_cat = dplyr::case_when(
        From != To & !grepl("^GAP", From) & !grepl("^GAP", To) &
          bayes_factor >= 10  & bayes_factor < 30  ~ "10–30 (strong)",
        From != To & !grepl("^GAP", From) & !grepl("^GAP", To) &
          bayes_factor >= 30  & bayes_factor < 100 ~ "30–100 (very strong)",
        From != To & !grepl("^GAP", From) & !grepl("^GAP", To) &
          bayes_factor >= 100                      ~ ">100 (decisive)",
        TRUE ~ NA_character_
      ),
      BF_star = dplyr::case_when(
        BF_cat == "10–30 (strong)"       ~ "*",
        BF_cat == "30–100 (very strong)" ~ "**",
        BF_cat == ">100 (decisive)"      ~ "***",
        TRUE ~ ""
      )
    )
  
  df <- df %>% filter(From %in% nodes, To %in% nodes)
  df$From <- factor(df$From, levels = nodes)
  df$To   <- factor(df$To,   levels = nodes)
  
  df_within <- df %>%
    filter(!grepl("^GAP", From), !grepl("^GAP", To)) %>%
    mutate(grp_from = get_group(From), grp_to = get_group(To)) %>%
    filter(grp_from == grp_to)
  
  p <- ggplot() +
    scale_x_discrete(drop = FALSE, labels = rename_labels, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, labels = rename_labels, expand = c(0, 0)) +
    labs(x = "Sink", y = "Source") +
    panel_theme +
    theme(
      axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 11, margin = margin(r = 2)),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold",
                                  margin = margin(r = 40))
    )
  
  for (i in seq_along(groups_keep)) {
    g      <- groups_keep[i]
    subdat <- df %>% filter(get_group(From) == g)
    
    p <- p +
      geom_tile(
        data = subdat,
        aes(x = To, y = From, fill = Rate_plot),
        color = "black", linewidth = 0.3
      ) +
      scale_fill_gradient(
        name     = if (i == 1) "Transition rate" else NULL,
        low      = group_cols[[g]][1],
        high     = group_cols[[g]][2],
        limits   = c(0, global_max),
        breaks   = c(0, 1, 2, 3),
        na.value = "white",
        guide    = guide_colorbar(
          order          = i,
          direction      = "horizontal",
          title.position = "top",
          barwidth       = LEG_BARWIDTH,
          barheight      = LEG_BARHEIGHT,
          ticks          = (i == length(groups_keep)),
          label          = (i == length(groups_keep)),
          label.position = "bottom"
        )
      )
    
    if (i != length(groups_keep)) {
      p <- p + ggnewscale::new_scale_fill()
    }
  }
  
  p <- p +
    geom_tile(
      data = df_within,
      aes(x = To, y = From),
      fill = NA, colour = "black", linewidth = 0.4
    ) +
    geom_text(
      data = df %>% filter(BF_star != ""),
      aes(x = To, y = From, label = BF_star),
      size = 4, fontface = "bold"
    )
  
  bf_legend_df <- tibble(
    To     = factor(levels(df$To)[1],   levels = levels(df$To)),
    From   = factor(levels(df$From)[1], levels = levels(df$From)),
    BF_cat = factor(
      c("10–30 (strong)", "30–100 (very strong)", ">100 (decisive)"),
      levels = c("10–30 (strong)", "30–100 (very strong)", ">100 (decisive)")
    )
  )
  
  p <- p +
    geom_text(
      data      = bf_legend_df,
      aes(x = To, y = From, colour = BF_cat, label = "a"),
      size = 6, fontface = "bold", alpha = 0, show.legend = TRUE
    ) +
    scale_colour_manual(
      name   = "Bayes factor",
      values = c("black", "black", "black"),
      guide  = guide_legend(
        order          = length(groups_keep) + 1,
        title.position = "top",
        direction      = "vertical",
        nrow           = 3,
        byrow          = TRUE,
        override.aes   = list(alpha = 1, label = c("*", "**", "***"))
      )
    )
  
  p <- p +
    annotate(
      "text",
      x          = -0.6,
      y          = side_label_df$ymid,
      label      = side_label_df$label,
      angle      = 90,
      hjust      = 0.5,
      vjust      = 1,
      size       = 4,
      fontface   = "bold",
      lineheight = 0.8
    ) +
    coord_cartesian(clip = "off")
  
  print(p)
  
  out_file <- paste0("Regions_all_HG_Rates_heatmap_median_", tag, ".png")
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 22,
    height   = 19,
    bg       = "white",
    dpi      = 300
  )
  message("✓ ", out_file, " saved.")
}

