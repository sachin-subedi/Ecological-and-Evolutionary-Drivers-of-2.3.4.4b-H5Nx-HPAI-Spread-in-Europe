options(scipen = 999)

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(purrr)


setwd("GLM")
base_dir <- "~/pipeline"

glm_dirs <- c(
  Subsample1 = file.path(base_dir, "subsampled1/combined/reproduce/GLM"),
  Subsample2 = file.path(base_dir, "subsampled2/combined/reproduce/GLM"),
  Subsample3 = file.path(base_dir, "subsampled3/combined/reproduce/GLM")
)

read_glm_subsample <- function(sub_name, glm_dir) {
  file_paths <- list(
    Set1a = file.path(glm_dir, "Set1a", sprintf("%s_Set1a_GLMSummary.txt", sub_name)),
    Set2  = file.path(glm_dir, "Set2",  sprintf("%s_Set2_GLMSummary.txt",  sub_name)),
    Set3  = file.path(glm_dir, "Set3",  sprintf("%s_Set3_GLMSummary.txt",  sub_name)),
    Set4  = file.path(glm_dir, "Set4",  sprintf("%s_Set4_GLMSummary.txt",  sub_name)),
    Set5  = file.path(glm_dir, "Set5",  sprintf("%s_Set5_GLMSummary.txt",  sub_name))
  )
  
  imap_dfr(
    file_paths,
    ~ read_tsv(.x, show_col_types = FALSE) %>%
      mutate(Set = .y,
             Subsample = sub_name)
  )
}

combined_glm <- map2_dfr(
  names(glm_dirs),
  glm_dirs,
  read_glm_subsample
)

combined_glm$Variable <- factor(
  combined_glm$Variable,
  levels = c(
    "Case Counts Origin",                 "Case Counts Destination",
    "Sample Size Origin",                 "Sample Size Destination",
    "Anatidae Counts Origin",             "Anatidae Counts Destination",
    "Accipitridae Counts Origin",         "Accipitridae Counts Destination",
    "Laridae Counts Origin",              "Laridae Counts Destination",
    "Live Poultry (NetWeight)",           "Live Poultry (Value)",
    "Road Transported Goods Origin",      "Road Transported Goods Destination",
    "Poultry Population Origin",          "Poultry Population Destination",
    "Animal Population Origin",           "Animal Population Destination",
    "Wetland Areas Origin",               "Wetland Areas Destination",
    "Forest Areas Origin",                "Forest Areas Destination",
    "Agricultural Land Origin",           "Agricultural Land Destination",
    "Water Bodies Areas Origin",          "Water Bodies Areas Destination",
    "Temperature Origin",                 "Temperature Destination",
    "Rainfall Origin",                    "Rainfall Destination",
    "Humidity Origin",                    "Humidity Destination",
    "Wind Speed Origin",                  "Wind Speed Destination",
    "High Vegetation Origin",             "High Vegetation Destination",
    "Low Vegetation Origin",              "Low Vegetation Destination",
    "Anatidae Migration"
  )
)

sub_levels <- c("Subsample1", "Subsample2", "Subsample3")
combined_glm$Subsample <- factor(combined_glm$Subsample, levels = sub_levels)

var_levels <- levels(combined_glm$Variable)
varsub_levels <- unlist(
  lapply(var_levels, function(v) paste(v, sub_levels, sep = " – "))
)

combined_glm <- combined_glm %>%
  mutate(
    VarSub = paste(Variable, Subsample, sep = " – ")
  )

combined_glm$VarSub <- factor(
  combined_glm$VarSub,
  levels = rev(varsub_levels)
)


combined_glm <- combined_glm %>%
  mutate(
    Significant = case_when(
      CoefLowerHPD > 0 ~ "Positive",
      CoefUpperHPD < 0 ~ "Negative",
      TRUE             ~ "Non-Significant"
    )
  )


sub_cols <- c(
  "Subsample1" = "#F39B7F",
  "Subsample2" = "#1B9E77",
  "Subsample3" = "#3C5488"
)


shape_map <- c(
  "Positive"        = 16,
  "Negative"        = 16,
  "Non-Significant" = 16
)

effect_size_plot <- ggplot(
  combined_glm,
  aes(x = CoefMedian, y = VarSub, colour = Subsample)
) +
  geom_point(aes(size = pp, shape = Significant)) +
  geom_errorbarh(
    aes(xmin = CoefLowerHPD, xmax = CoefUpperHPD),
    height = 0.2,
    colour = "black"
  ) +
  geom_vline(xintercept = 0, colour = "black") +
  scale_colour_manual(values = sub_cols, name = "Subsample", guide = "none") +
  scale_shape_manual(values = shape_map, guide = "none") +
  scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
  scale_y_discrete(labels = function(x) {
    vars <- gsub(" – Subsample[123]", "", x)
    subs <- sub(".* – (Subsample[123])$", "\\1", x)
    ifelse(subs == "Subsample2", vars, "")
  }) +
  labs(x = "Conditional effect size", y = NULL) +
  theme_classic(base_family = "Times New Roman") +
  theme(
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black"),
    axis.text.y  = element_text(size = 16, face = "bold", colour = "black"),
    axis.line.x  = element_line(colour = "black"),
    axis.line.y  = element_line(colour = "black"),
    axis.ticks   = element_line(colour = "black")
  )

prior  <- 0.05
bf_ref <- c(10, 30, 100)
pp_ref <- (bf_ref * prior) / (bf_ref * prior + 1 - prior)

posterior_prob_plot <- ggplot(
  combined_glm,
  aes(x = pp, y = VarSub, fill = Subsample)
) +
  geom_col(width = 0.9, colour = "black") +   
  geom_vline(xintercept = pp_ref, linetype = "dashed", colour = "darkblue") +
  annotate("text", x = pp_ref[1], y = Inf, label = "BF = 10–30",
           vjust = -0.6, size = 4.5, family = "Times New Roman") +
  annotate("text", x = pp_ref[2], y = Inf, label = "BF = 30–100",
           vjust = -0.6, size = 4.5, family = "Times New Roman") +
  annotate("text", x = pp_ref[3], y = Inf, label = "BF > 100",
           vjust = -0.6, hjust = 1, size = 4.5, family = "Times New Roman") +
  scale_fill_manual(values = sub_cols, name = "Subsample") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = "Posterior probability", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Times New Roman") +
  theme(
    axis.title.x   = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x    = element_text(size = 12, colour = "black"),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    axis.line.x    = element_line(colour = "black"),
    axis.line.y    = element_blank(),
    plot.margin     = margin(t = 24, r = 32)
  )


combined_plot <- effect_size_plot + posterior_prob_plot +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(
    plot.margin      = margin(10, 20, 10, 20),
    axis.text.x      = element_text(colour = "black"),
    legend.position  = "bottom",                    
    legend.direction = "horizontal",
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 12),
    legend.key.width = unit(1.5, "lines")
  )

combined_plot

ggsave(
  filename = "AllSubsamples_GeoCluster_GLM_cluster_summary_plot_NS.png",
  plot     = combined_plot,
  width    = 20,
  height   = 18,
  dpi      = 300,
  units    = "in"
)

glm_summary <- combined_glm %>%
  dplyr::select(Subsample, Variable, CoefMedian, CoefLowerHPD, CoefUpperHPD, pp, BF) %>%
  arrange(desc(BF))

write_csv(glm_summary,
          "AllSubsamples_GeoCluster_GLM_glm_summary_sortedBF.csv")

