suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
  library(lubridate)
  library(scales)
})

csv_path <- "filtered_data_combined_Anatidae.csv"

sumw_path <- "/NS_biased/season_migratory_raw/Anatidae_migratory_SUMW.tsv"

out_dir <- "NS_biased/season_migratory_raw/ebird_density_check/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

years    <- 2016:2025
DIAG_EPS <- 1e-11

cluster_region_lookup <- tribble(
  ~Country,                    ~Region,
  "Albania",                   "HC_Cluster_2_Mediterranean",
  "Austria",                   "HC_Cluster_2_Alpine",
  "Belgium",                   "HC_Cluster_1_Atlantic",
  "Belarus",                   "HC_Cluster_2_Alpine",
  "Bosnia and Herz.",          "HC_Cluster_2_Alpine",
  "Bulgaria",                  "HC_Cluster_2_Continental",
  "Croatia",                   "HC_Cluster_2_Continental",
  "Cyprus",                    "HC_Cluster_2_Mediterranean",
  "Czechia",                   "HC_Cluster_2_Continental",
  "Denmark",                   "HC_Cluster_2_Continental",
  "Estonia",                   "HC_Cluster_3_Boreal",
  "Finland",                   "HC_Cluster_3_Boreal",
  "France",                    "HC_Cluster_1_Atlantic",
  "Germany",                   "HC_Cluster_1_Continental",
  "Greece",                    "HC_Cluster_2_Mediterranean",
  "Hungary",                   "HC_Cluster_2_Pannonian",
  "Iceland",                   "HC_Cluster_1_Atlantic",
  "Ireland",                   "HC_Cluster_1_Atlantic",
  "Italy",                     "HC_Cluster_2_Mediterranean",
  "Kosovo",                    "HC_Cluster_2_Continental",
  "Latvia",                    "HC_Cluster_3_Boreal",
  "Lithuania",                 "HC_Cluster_3_Boreal",
  "Luxembourg",                "HC_Cluster_1_Continental",
  "Moldova",                   "HC_Cluster_2_Continental",
  "Montenegro",                "HC_Cluster_2_Alpine",
  "Netherlands",               "HC_Cluster_1_Atlantic",
  "North Macedonia",           "HC_Cluster_2_Continental",
  "Norway",                    "HC_Cluster_3_Alpine",
  "Poland",                    "HC_Cluster_2_Continental",
  "Portugal",                  "HC_Cluster_4_Mediterranean",
  "Romania",                   "HC_Cluster_2_Continental",
  "Serbia",                    "HC_Cluster_2_Continental",
  "Slovakia",                  "HC_Cluster_2_Alpine",
  "Slovenia",                  "HC_Cluster_2_Continental",
  "Spain",                     "HC_Cluster_4_Mediterranean",
  "Sweden",                    "HC_Cluster_3_Boreal",
  "Switzerland",               "HC_Cluster_1_Alpine",
  "Ukraine",                   "HC_Cluster_2_Continental",
  "United Kingdom",            "HC_Cluster_1_Atlantic"
)

clusters <- c(
  "HC_Cluster_1_Alpine", "HC_Cluster_1_Atlantic", "HC_Cluster_1_Continental",
  "HC_Cluster_2_Alpine", "HC_Cluster_2_Continental", "HC_Cluster_2_Mediterranean",
  "HC_Cluster_2_Pannonian", "HC_Cluster_3_Alpine", "HC_Cluster_3_Boreal",
  "HC_Cluster_4_Mediterranean"
)

cluster_pal <- c(
  "HC_Cluster_1_Alpine"        = "#CC79A7",
  "HC_Cluster_1_Atlantic"      = "#CD5C5C",
  "HC_Cluster_1_Continental"   = "#BCBD22",
  "HC_Cluster_2_Alpine"        = "#26A69A",
  "HC_Cluster_2_Continental"   = "#0072B2",
  "HC_Cluster_2_Mediterranean" = "lightslategray",
  "HC_Cluster_2_Pannonian"     = "#6B4C1B",
  "HC_Cluster_3_Alpine"        = "#7570B3",
  "HC_Cluster_3_Boreal"        = "#006D2C",
  "HC_Cluster_4_Mediterranean" = "#56B4E9"
)

cluster_labels <- c(
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

lab_gc <- function(x) cluster_labels[x]

season_from_doy <- function(d) {
  d <- as.integer(d)
  ifelse(d >= 46  & d < 135, "spring",
         ifelse(d >= 258 & d < 319, "fall", "nonmig"))
}

message("• Loading eBird data...")

birds <- read_csv(csv_path,
                  col_types = cols(
                    startDayOfYear   = col_double(),
                    individualCount  = col_double(),
                    decimalLatitude  = col_double(),
                    decimalLongitude = col_double(),
                    genus            = col_character(),
                    .default         = col_guess()
                  ), show_col_types = FALSE) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = 4326, remove = FALSE)

europe <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(name != "Russia") %>%
  rename(Country = name) %>%
  st_make_valid()

birds <- st_join(birds, europe["Country"], join = st_within, left = FALSE) %>%
  left_join(cluster_region_lookup, by = "Country") %>%
  filter(!is.na(Region))

cols_b <- names(birds)
if (all(c("year", "month", "day") %in% cols_b)) {
  birds$year  <- as.integer(birds$year)
  birds$month <- pmax(1L, pmin(12L, as.integer(birds$month)))
  d_raw       <- suppressWarnings(as.integer(birds$day))
  birds$day   <- ifelse(is.na(d_raw) | d_raw < 1L | d_raw > 31L, 1L, d_raw)
  birds$date  <- make_date(year = birds$year, month = birds$month, day = birds$day)
} else {
  birds$year           <- as.integer(birds$year)
  birds$startDayOfYear <- pmax(1, pmin(366, as.integer(round(birds$startDayOfYear))))
  birds$date <- ymd(paste0(birds$year, "-01-01")) + days(birds$startDayOfYear - 1)
}

birds$doy    <- yday(birds$date)
birds$season <- season_from_doy(birds$doy)
birds        <- birds %>% filter(year %in% years)

message(sprintf("  ✓ %s records loaded across %s countries",
                format(nrow(birds), big.mark = ","),
                n_distinct(birds$Country)))

message("• Summarising observation density...")

obs_density <- birds %>%
  st_drop_geometry() %>%
  filter(season %in% c("spring", "fall")) %>%
  group_by(Region, year) %>%
  summarise(
    n_records = n(),
    n_cells   = n_distinct(paste(round(decimalLongitude, 1),
                                 round(decimalLatitude,  1))),
    n_genera  = n_distinct(genus),
    .groups   = "drop"
  ) %>%
  mutate(Region = factor(Region, levels = clusters))

p1 <- ggplot(obs_density,
             aes(x = factor(year), y = n_records, fill = Region)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = cluster_pal,
    name   = NULL,
    labels = lab_gc
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    x     = "Year",
    y     = "eBird Anatidae records (Spring + Fall)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    legend.text     = element_text(size = 10),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.title      = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 3, title = NULL))

ggsave(file.path(out_dir, "Fig_S_ebird_density_by_year.png"),
       p1, width = 12, height = 6, dpi = 300, bg = "white")
message("  ✓ Fig_S_ebird_density_by_year.png")

p2 <- ggplot(obs_density,
             aes(x = Region, y = n_records, fill = Region)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "grey20") +
  scale_fill_manual(values = cluster_pal, guide = "none") +
  scale_x_discrete(labels = lab_gc) +
  scale_y_continuous(labels = comma) +
  coord_flip() +
  labs(
    x = NULL,
    y = "eBird Anatidae records per year (Spring + Fall)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Fig_S_ebird_density_boxplot.png"),
       p2, width = 9, height = 6, dpi = 300, bg = "white")
message("  ✓ Fig_S_ebird_density_boxplot.png")

obs_country <- birds %>%
  st_drop_geometry() %>%
  filter(season %in% c("spring", "fall"), year %in% years) %>%
  group_by(Country) %>%
  summarise(mean_records = n() / length(years), .groups = "drop")

country_map <- europe %>%
  left_join(obs_country, by = "Country")

p3 <- ggplot(country_map) +
  geom_sf(aes(fill = mean_records), color = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    name      = "Mean annual\nrecords\n(Spring + Fall)",
    labels    = comma,
    na.value  = "grey85"
  ) +
  coord_sf(xlim = c(-25, 41), ylim = c(34, 81), expand = FALSE) +
  labs(title = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    panel.grid.major = element_line(color = "grey92")
  )

ggsave(file.path(out_dir, "Fig_S_ebird_density_map.png"),
       p3, width = 8, height = 7, dpi = 300, bg = "white")
message("  ✓ Fig_S_ebird_density_map.png")

message("• Running correlation: eBird density vs SUMW...")

M_SUM <- read_tsv(sumw_path, show_col_types = FALSE) %>%
  column_to_rownames("A") %>%
  as.matrix()

sumw_long <- as_tibble(M_SUM, rownames = "From") %>%
  pivot_longer(-From, names_to = "To", values_to = "SUMW") %>%
  filter(From != To, SUMW > DIAG_EPS)

obs_total <- obs_density %>%
  group_by(Region) %>%
  summarise(total_records = sum(n_records), .groups = "drop") %>%
  rename(From = Region)

corr_df <- sumw_long %>%
  left_join(obs_total, by = "From") %>%
  filter(!is.na(total_records))

cor_test <- cor.test(corr_df$total_records, corr_df$SUMW, method = "spearman")

cat(sprintf(
  "\n══════════════════════════════════════\n  Spearman rho = %.3f\n  p-value      = %.4f\n  n pairs      = %d\n══════════════════════════════════════\n",
  cor_test$estimate, cor_test$p.value, nrow(corr_df)
))

writeLines(
  sprintf(
    "eBird density vs migration network strength (SUMW)\nSpearman rho = %.3f\np-value      = %.4f\nn pairs      = %d",
    cor_test$estimate, cor_test$p.value, nrow(corr_df)
  ),
  file.path(out_dir, "correlation_result.txt")
)

p4 <- ggplot(corr_df,
             aes(x = total_records, y = SUMW, color = From)) +
  geom_smooth(
    method      = "lm",
    se          = TRUE,
    color       = "grey30",
    linewidth   = 0.8,
    linetype    = "dashed",
    alpha       = 0.12,
    inherit.aes = FALSE,
    aes(x = total_records, y = SUMW)
  ) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    values = cluster_pal,
    name   = NULL,
    labels = lab_gc
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  annotate(
    "text",
    x      = Inf, y = Inf,
    label  = sprintf("Spearman \u03c1 = %.2f\np = %.3f",
                     cor_test$estimate, cor_test$p.value),
    hjust  = 1.1, vjust = 1.5,
    size   = 5, fontface = "italic"
  ) +
  labs(
    x = "Total eBird Anatidae records — source Region (Spring + Fall, 2016–2025)",
    y = "Mean migration corridor strength (SUMW)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    legend.text     = element_text(size = 10),
    axis.title      = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 3, title = NULL))

ggsave(file.path(out_dir, "Fig_S_ebird_vs_network.png"),
       p4, width = 9, height = 7, dpi = 300, bg = "white")

