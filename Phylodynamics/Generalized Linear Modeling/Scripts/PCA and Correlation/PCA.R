library(tidyverse)
library(FactoMineR)
library(factoextra)
library(plotly)
library(corrplot)
library(car)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/GLM_code_predictors/PCA/")
predictor_dir <- "~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/GLM_code_predictors/PCA/"
tsv_files     <- list.files(predictor_dir, pattern = "\\.tsv$", full.names = TRUE)

data_list <- lapply(tsv_files, function(f) read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
names(data_list) <- basename(tsv_files) %>% str_remove("\\.tsv$")

key_col <- colnames(data_list[[1]])[1]

merged <- reduce(data_list, left_join, by = key_col)

key_col <- colnames(merged)[1]
df_num  <- merged %>%
  select(-all_of(key_col)) %>%
  mutate(across(where(is.character), readr::parse_number)) %>%
  select(where(is.numeric))

colnames(df_num) <- c(
  "Accipitridae Counts",      # accipitridae_counts
  "Agricultural Land Areas",  # AgriculturalLand_areas
  "Anatidae Counts",          # anatidae_counts
  "Forest Areas",             # ForestAreas_areas
  "Animal Population",        # Animal_Population
  "Case Counts",              # case_counts
  "Humidity",                 # humidity
  "High Vegetation",          # LAI_HV
  "Low Vegetation",           # LAI_LV
  "Livestock Density Index",  # livestock_density_index
  "Poultry Population",       # poultry_population
  "Rainfall",                 # rainfall
  "Road Transported Goods",   # road_transported_goods
  "Temperature",              # temperature
  "Wind Speed",               # windspeed
  "Laridae Counts",           # laridae_counts
  "Water Bodies Areas",       # WaterBodies_areas
  "Wetland Areas"             # Wetlands_areas
)

na_frac <- colMeans(is.na(df_num))
df_num  <- df_num[, na_frac < 0.5, drop = FALSE]
nzv     <- sapply(df_num, function(x) sd(x, na.rm = TRUE) > 0)
df_num  <- df_num[, nzv, drop = FALSE]

df_scaled <- scale(df_num)
pca <- prcomp(df_scaled, center = FALSE, scale. = FALSE)

fviz_eig(pca, addlabels = TRUE, barfill = "grey70", barcolor = "grey40",
         linecolor = "grey20", ggtheme = theme_minimal()) +
  labs(title = "PCA: Variance Explained")

fviz_pca_biplot(
  pca,
  label = "var",         # show variable loadings labels
  geom.ind = "point",    # points for GeoClusters (observations)
  col.ind = "grey35",
  repel = TRUE
) + theme_minimal() + labs(title = "PCA Biplot: PC1 vs PC2")

loadings <- as.data.frame(pca$rotation[, 1:3])
loadings$Variable <- rownames(loadings)

plot_ly(
  loadings,
  x = ~PC1, y = ~PC2, z = ~PC3,
  text = ~Variable, type = "scatter3d", mode = "text+markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "PCA Loadings: PC1 vs PC2 vs PC3",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

loadings <- as.data.frame(pca$rotation[, 1:3])
loadings$Variable <- rownames(loadings)

plot_ly(
  loadings,
  x = ~PC1, y = ~PC2, z = ~PC3,
  text = ~Variable,
  type = "scatter3d",
  mode = "markers+text",
  textposition = "top center",
  marker = list(
    size = 6,
    line = list(width = 1, color = "black")
  )
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = "PC1",
        showbackground = TRUE,
        backgroundcolor = "rgb(240,240,240)",
        gridcolor       = "rgb(210,210,210)",
        zerolinecolor   = "rgb(180,180,180)"
      ),
      yaxis = list(
        title = "PC2",
        showbackground = TRUE,
        backgroundcolor = "rgb(240,240,240)",
        gridcolor       = "rgb(210,210,210)",
        zerolinecolor   = "rgb(180,180,180)"
      ),
      zaxis = list(
        title = "PC3",
        showbackground = TRUE,
        backgroundcolor = "rgb(240,240,240)",
        gridcolor       = "rgb(210,210,210)",
        zerolinecolor   = "rgb(180,180,180)"
      ),
      bgcolor    = "rgb(255,255,255)",
      aspectmode = "cube",       
      camera     = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
    )
  )

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)


loadings <- as.data.frame(pca$rotation[, 1:3]) %>%
  rownames_to_column("Variable")   # PC1, PC2, PC3 + Variable


project_3d <- function(x, y, z, theta = 45, phi = 25) {
  th <- theta * pi / 180
  ph <- phi   * pi / 180
  
  # rotate around Z
  x1 <-  x * cos(th) + y * sin(th)
  y1 <- -x * sin(th) + y * cos(th)
  z1 <-  z
  
  # tilt around X
  x2 <- x1
  y2 <- y1 * cos(ph) + z1 * sin(ph)
  
  data.frame(x_proj = x2, y_proj = y2)
}

theta_view <- 45
phi_view   <- 25

proj        <- project_3d(loadings$PC1, loadings$PC2, loadings$PC3,
                          theta = theta_view, phi = phi_view)
loadings_2d <- bind_cols(loadings, proj)

L <- max(abs(unlist(loadings[, c("PC1", "PC2", "PC3")]))) * 1.3

axes_raw <- tibble(
  axis = rep(c("PC1", "PC2", "PC3"), each = 2),
  x    = c(0, L, 0, 0, 0, 0),
  y    = c(0, 0, 0, L, 0, 0),
  z    = c(0, 0, 0, 0, 0, L)
)

axes_proj <- bind_cols(
  axes_raw,
  project_3d(axes_raw$x, axes_raw$y, axes_raw$z,
             theta = theta_view, phi = phi_view)
)

axes_seg <- axes_proj %>%
  group_by(axis) %>%
  summarise(
    x    = first(x_proj),
    y    = first(y_proj),
    xend = last(x_proj),
    yend = last(y_proj),
    .groups = "drop"
  )

p_3d_gg <- ggplot() +
  geom_segment(
    data = axes_seg,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.6,
    colour = "grey40",
    arrow = arrow(length = unit(0.18, "cm"))
  ) +
  geom_text(
    data = axes_seg,
    aes(x = xend, y = yend, label = axis),
    hjust = -0.1, vjust = -0.1,
    size = 3.2
  ) +
  geom_point(
    data = loadings_2d,
    aes(x = x_proj, y = y_proj),
    size = 2.8
  ) +
  geom_text_repel(
    data = loadings_2d,
    aes(x = x_proj, y = y_proj, label = Variable),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    segment.colour = "grey70"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "grey96", colour = NA),
    plot.background  = element_rect(fill = "white",   colour = NA),
    panel.grid       = element_line(colour = "grey88"),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank()
  )

p_3d_gg

ggsave(
  "PCA_loadings_3D_style_ggplot.png",
  plot = p_3d_gg,
  width = 10,
  height = 8,
  bg = "transparent",
  dpi = 300
)

cor_mat <- cor(df_num, use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", type = "upper", tl.col = "black",
         tl.cex = 1, addCoef.col = NA, mar = c(0,0,1,0),
         title = "Predictor Correlations")


library(corrplot)

cor_mat <- cor(df_num, use = "pairwise.complete.obs")
png("Predictor_Correlations.png",
   width = 2800, height = 1600, res = 300)

corrplot(
  cor_mat,
  method = "color",
  type   = "upper",
  tl.col = "black",
  tl.cex = 0.7,
  tl.srt = 45,
  addCoef.col = "black",
  number.cex = 0.5,
  number.digits = 2,
  mar = c(0, 0, 1, 0)
)

dev.off()


df_num_numeric <- df_num[, vapply(df_num, is.numeric, logical(1)), drop = FALSE]

desc_stats <- purrr::imap_dfr(df_num_numeric, ~ {
  x <- .x
  tibble::tibble(
    Variable = .y,
    Mean = mean(x, na.rm = TRUE),
    SD   = stats::sd(x, na.rm = TRUE),
    Min  = min(x, na.rm = TRUE),
    Max  = max(x, na.rm = TRUE)
  )
}) %>%
  dplyr::mutate(
    dplyr::across(c(Mean, SD, Min, Max), ~ round(.x, 3)),
    Range = paste0(Min, "–", Max)
  ) %>%
  dplyr::relocate(Variable, Mean, SD, Range)

# Check in console
print(desc_stats, n = Inf)

# Export to CSV
readr::write_csv(desc_stats,
                 "Predictor_descriptive_stats_mean_SD_range.csv")

