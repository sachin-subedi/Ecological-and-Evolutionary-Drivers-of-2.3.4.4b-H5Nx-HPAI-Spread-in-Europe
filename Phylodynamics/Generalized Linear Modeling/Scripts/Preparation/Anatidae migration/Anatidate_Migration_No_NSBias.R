library(tidyverse)
library(sf)
library(rnaturalearth)
library(FNN)
library(geodist)
library(lubridate)
library(igraph)
library(readr)
library(scales)
library(grid)

sf_use_s2(TRUE)

csv_path <- "filtered_data_combined_Anatidae.csv"
out_dir  <- "unbiased"
years    <- 2016:2025

DIAG_EPS <- 1e-11

cluster_region_lookup <- tribble(
  ~Country,                  ~GeoCluster,
  "Albania","GeoCluster_One","Austria","GeoCluster_Five","Belgium","GeoCluster_Three",
  "Bosnia and Herzegovina","GeoCluster_Five","Bulgaria","GeoCluster_One","Croatia","GeoCluster_Five",
  "Cyprus","GeoCluster_One","Czechia","GeoCluster_Five","Denmark","GeoCluster_Three",
  "Estonia","GeoCluster_Four","Finland","GeoCluster_Four","France","GeoCluster_Three",
  "Germany","GeoCluster_Three","Greece","GeoCluster_One","Hungary","GeoCluster_Two",
  "Iceland","GeoCluster_Three","Ireland","GeoCluster_Three","Italy","GeoCluster_Five",
  "Kosovo","GeoCluster_One","Latvia","GeoCluster_Four","Lithuania","GeoCluster_Four",
  "Luxembourg","GeoCluster_Three","Moldova","GeoCluster_Two","Netherlands","GeoCluster_Three",
  "North Macedonia","GeoCluster_One","Norway","GeoCluster_Four","Poland","GeoCluster_Two",
  "Portugal","GeoCluster_Three","Romania","GeoCluster_One","Serbia","GeoCluster_One",
  "Slovakia","GeoCluster_Two","Slovenia","GeoCluster_Five","Spain","GeoCluster_Three",
  "Sweden","GeoCluster_Four","Switzerland","GeoCluster_Three","Ukraine","GeoCluster_Two",
  "United Kingdom","GeoCluster_Three"
)

clusters <- c("GeoCluster_One","GeoCluster_Two","GeoCluster_Three","GeoCluster_Four","GeoCluster_Five")

cluster_pretty <- c(
  "GeoCluster_One"   = "South-Eastern",
  "GeoCluster_Two"   = "Central-Eastern",
  "GeoCluster_Three" = "Central-Western",
  "GeoCluster_Four"  = "Northern",
  "GeoCluster_Five"  = "Central-Southern"
)

lab_gc <- function(x){
  x <- as.character(x)
  out <- unname(cluster_pretty[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

zone_centers <- tibble::tibble(
  GeoCluster = factor(clusters, levels = clusters),
  x = c(25.7, 23.0, -1.5, 10.8, 12.5),
  y = c(45.9, 49.5, 52.5, 60.5, 45.5)
) %>% mutate(Region = lab_gc(GeoCluster))

season_from_doy <- function(d) {
  d <- as.integer(d)
  ifelse(d >= 46  & d < 135, "spring",
         ifelse(d >= 258 & d < 319, "fall", "nonmig"))
}

ensure_dir <- function(p){ if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE) }

to_square <- function(df_AB_val, diag_eps = DIAG_EPS){
  grid <- tidyr::expand_grid(A = clusters, B = clusters)
  mdf  <- grid %>% left_join(df_AB_val, by = c("A","B")) %>%
    mutate(val = replace_na(val, 0)) %>%
    mutate(A = factor(A, levels = clusters),
           B = factor(B, levels = clusters)) %>%
    arrange(A, B)
  wide <- mdf %>% select(A,B,val) %>% pivot_wider(names_from = B, values_from = val) %>% arrange(A)
  mm <- as.matrix(wide[,-1, drop = FALSE]); rownames(mm) <- wide$A
  diag(mm) <- diag_eps
  mm
}

write_matrix_tsv <- function(mat, file_path){
  out <- tibble(A = rownames(mat)) %>% bind_cols(as_tibble(mat, .name_repair = "minimal"))
  readr::write_tsv(out, file_path)
}

jacc_fast <- function(li, lj){
  if (is.null(li) || is.null(lj)) return(0)
  inter <- length(intersect(li, lj)); uni <- length(union(li, lj))
  if (uni == 0) 0 else inter / uni
}

message("• Loading Anatidae points …")
birds <- readr::read_csv(
  csv_path,
  col_types = cols(
    startDayOfYear   = col_double(),
    individualCount  = col_double(),
    decimalLatitude  = col_double(),
    decimalLongitude = col_double(),
    family           = col_character(),
    genus            = col_character(),
    .default         = col_guess()
  ),
  show_col_types = FALSE
) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(individualCount)) %>%
  st_as_sf(coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE)

europe <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(name != "Russia") %>%
  rename(Country = name) %>%
  st_make_valid()

birds <- st_join(birds, europe["Country"], join = st_within, left = FALSE) %>%
  left_join(cluster_region_lookup, by = "Country") %>%
  filter(!is.na(GeoCluster))

cols <- names(birds)
if (all(c("year","month","day") %in% cols)) {
  birds$year  <- as.integer(birds$year)
  birds$month <- pmax(1L, pmin(12L, as.integer(birds$month)))
  d_raw <- suppressWarnings(as.integer(birds$day))
  birds$day  <- ifelse(is.na(d_raw) | d_raw < 1L | d_raw > 31L, 1L, d_raw)
  birds$date <- make_date(year = birds$year, month = birds$month, day = birds$day)
} else if (all(c("year","startDayOfYear") %in% cols)) {
  birds$year <- as.integer(birds$year)
  birds$startDayOfYear <- pmax(1, pmin(366, as.integer(round(birds$startDayOfYear))))
  birds$date <- ymd(paste0(birds$year, "-01-01")) + days(birds$startDayOfYear - 1)
} else {
  stop("No usable date fields. Need either (year,month,day) or (year,startDayOfYear).")
}

birds$doy    <- yday(birds$date)
birds$season <- season_from_doy(birds$doy)
birds        <- birds %>% filter(year %in% years)

eu_3035       <- st_transform(europe, 3035)
eu_union_3035 <- st_union(eu_3035)
europe_union  <- st_transform(eu_union_3035, 4326)

grid_05 <- st_make_grid(europe_union, cellsize = c(0.5,0.5), what = "polygons", square = TRUE) %>%
  st_sf(crs = 4326) %>%
  st_filter(europe_union, .predicate = st_intersects)

cc <- st_coordinates(st_centroid(st_geometry(grid_05)))
grid_05$grid_id <- paste0("N_", round(cc[, "Y"], 3), "_", round(cc[, "X"], 3))

birds_grid <- st_join(birds, grid_05["grid_id"], join = st_within, left = FALSE)

nodes_tbl <- birds_grid %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    lat = mean(decimalLatitude), lon = mean(decimalLongitude),
    Country = names(which.max(table(Country))),
    GeoCluster = names(which.max(table(GeoCluster))),
    .groups = "drop"
  )
stopifnot(nrow(nodes_tbl) >= 2)

nodes_sf <- merge(grid_05["grid_id"], nodes_tbl, by = "grid_id", all = FALSE)

set.seed(42)
coords_mat <- nodes_tbl %>% select(lon, lat) %>% as.matrix()
n_nodes <- nrow(coords_mat)
k_use   <- max(1, min(20, n_nodes - 1))   # k up to 20
knn     <- FNN::get.knn(coords_mat, k = k_use)

edges_idx <- tibble(
  src = rep(seq_len(n_nodes), each = ncol(knn$nn.index)),
  dst = as.integer(as.vector(knn$nn.index))
) %>%
  filter(src != dst) %>%
  mutate(a = pmin(src, dst), b = pmax(src, dst)) %>%
  distinct(a, b) %>%
  transmute(src = a, dst = b)

lat_i <- nodes_tbl$lat[edges_idx$src]; lon_i <- nodes_tbl$lon[edges_idx$src]
lat_j <- nodes_tbl$lat[edges_idx$dst]; lon_j <- nodes_tbl$lon[edges_idx$dst]

DIST_ij  <- geodist::geodist(
  data.frame(x = lon_i, y = lat_i),
  data.frame(x = lon_j, y = lat_j),
  measure = "geodesic", paired = TRUE
) / 1000

A_ids <- nodes_tbl$grid_id[edges_idx$src]
B_ids <- nodes_tbl$grid_id[edges_idx$dst]

node_genus_by_year_season <- function(y, s){
  bys <- birds_grid %>% filter(year == y, season == s) %>% st_drop_geometry()
  if (nrow(bys) == 0) return(NULL)
  bys %>%
    group_by(grid_id) %>%
    summarise(
      genus_list = list(sort(unique(genus))),
      abund      = sum(individualCount, na.rm = TRUE),
      lat        = mean(decimalLatitude),
      GeoCluster = names(which.max(table(GeoCluster))),
      .groups = "drop"
    )
}

build_migratory_for_year <- function(y){
  eps <- 1e-6
  
  ng_s <- node_genus_by_year_season(y, "spring")
  spring_df <- NULL
  if (!is.null(ng_s)) {
    gl <- setNames(ng_s$genus_list, ng_s$grid_id)
    cl <- setNames(ng_s$GeoCluster, ng_s$grid_id)
    lt <- setNames(ng_s$lat,        ng_s$grid_id)
    
    keep <- (A_ids %in% ng_s$grid_id) & (B_ids %in% ng_s$grid_id)
    if (any(keep)) {
      src <- A_ids[keep]; dst <- B_ids[keep]
      lat_i_y <- lt[src];  lat_j_y <- lt[dst]
      Gs <- cl[src];       Gt <- cl[dst]
      
      J <- mapply(\(i,j) jacc_fast(gl[[i]], gl[[j]]), src, dst, USE.NAMES = FALSE)
      D <- as.numeric(DIST_ij[keep])
      W <- J / pmax(D, eps)
      
      A <- ifelse(lat_i_y <= lat_j_y, Gs, Gt)
      B <- ifelse(lat_i_y <= lat_j_y, Gt, Gs)
      
      spring_df <- tibble(A=A, B=B, W=W, J=J)
    }
  }
  
  ng_f <- node_genus_by_year_season(y, "fall")
  fall_df <- NULL
  if (!is.null(ng_f)) {
    gl <- setNames(ng_f$genus_list, ng_f$grid_id)
    cl <- setNames(ng_f$GeoCluster, ng_f$grid_id)
    lt <- setNames(ng_f$lat,        ng_f$grid_id)
    
    keep <- (A_ids %in% ng_f$grid_id) & (B_ids %in% ng_f$grid_id)
    if (any(keep)) {
      src <- A_ids[keep]; dst <- B_ids[keep]
      lat_i_y <- lt[src];  lat_j_y <- lt[dst]
      Gs <- cl[src];       Gt <- cl[dst]
      
      J <- mapply(\(i,j) jacc_fast(gl[[i]], gl[[j]]), src, dst, USE.NAMES = FALSE)
      D <- as.numeric(DIST_ij[keep])
      W <- J / pmax(D, eps)
      
      A <- ifelse(lat_i_y >= lat_j_y, Gs, Gt)
      B <- ifelse(lat_i_y >= lat_j_y, Gt, Gs)
      
      fall_df <- tibble(A=A, B=B, W=W, J=J)
    }
  }
  
  mig <- bind_rows(spring_df, fall_df)
  if (is.null(mig) || nrow(mig) == 0) return(NULL)
  
  mig %>%
    filter(A != B, J > 0, is.finite(W), W > 0) %>%
    group_by(A,B) %>%
    summarise(
      EDGE_COUNT = n(),
      SUMW       = sum(W, na.rm = TRUE),
      MEANW      = mean(W, na.rm = TRUE),
      .groups = "drop"
    ) %>% mutate(year = y, season = "migratory")
}

message("• Building per-year migratory summaries …")
per_year  <- purrr::compact(lapply(years, build_migratory_for_year))
pair_year <- bind_rows(per_year)
stopifnot(nrow(pair_year) > 0)

avg_by_season <- function(s, col){
  grid_yab <- tidyr::expand_grid(year = years, A = clusters, B = clusters) %>% filter(A != B)
  dat <- pair_year %>% filter(season == s) %>%
    right_join(grid_yab, by = c("year","A","B")) %>%
    mutate(
      EDGE_COUNT = replace_na(EDGE_COUNT, 0),
      SUMW       = replace_na(SUMW, 0),
      MEANW      = replace_na(MEANW, 0)
    )
  dat %>% group_by(A,B) %>% summarise(val = mean(.data[[col]], na.rm = TRUE), .groups = "drop")
}

ensure_dir(out_dir)
mig_dir <- file.path(out_dir, "season_migratory_raw_noNSbias")
ensure_dir(mig_dir)

den_avg <- avg_by_season("migratory", "EDGE_COUNT")
sum_avg <- avg_by_season("migratory", "SUMW")
men_avg <- avg_by_season("migratory", "MEANW")

M_DEN <- to_square(den_avg, diag_eps = DIAG_EPS)
M_SUM <- to_square(sum_avg, diag_eps = DIAG_EPS)
M_MEN <- to_square(men_avg, diag_eps = DIAG_EPS)

write_matrix_tsv(M_DEN, file.path(mig_dir, "Anatidae_migratory_EDGE_COUNT.tsv"))
write_matrix_tsv(M_SUM, file.path(mig_dir, "Anatidae_migratory_SUMW.tsv"))
write_matrix_tsv(M_MEN, file.path(mig_dir, "Anatidae_migratory_MEANW.tsv"))

message("✓ Wrote MIGRATORY matrices (RAW, no N–S bias) to: ", mig_dir)

sum_path <- file.path(mig_dir, "Anatidae_migratory_SUMW.tsv")
M_SUM0 <- readr::read_tsv(sum_path, show_col_types = FALSE) |>
  column_to_rownames("A") |>
  as.matrix()

stopifnot(all(rownames(M_SUM0) %in% clusters), all(colnames(M_SUM0) %in% clusters))
M_SUM0 <- M_SUM0[clusters, clusters, drop = FALSE]

hm_df <- as_tibble(M_SUM0, rownames = "A") |>
  pivot_longer(-A, names_to = "B", values_to = "val") |>
  mutate(
    A = factor(A, levels = clusters),
    B = factor(B, levels = clusters),
    val_plot = if_else(as.character(A) == as.character(B), NA_real_, val)
  )

v_pos <- hm_df$val_plot[is.finite(hm_df$val_plot) & hm_df$val_plot > 0]
low  <- min(v_pos, na.rm = TRUE)
high <- max(v_pos, na.rm = TRUE)
mid  <- exp(mean(log(v_pos)))      # geometric mean
brks <- c(low, mid, high)

fmt3 <- function(x){
  top <- max(v_pos, na.rm = TRUE)
  acc <- if (top < 0.1) 0.001 else if (top < 1) 0.01 else 0.1
  label_number(accuracy = acc)(x)
}

p_hm <- ggplot(hm_df, aes(B, A, fill = val_plot)) +
  geom_tile(color = "black", linewidth = 0.6) +
  scale_fill_gradientn(
    colours = c("#e0f3f8", "#abd9e9", "#2166ac"),
    trans   = pseudo_log_trans(sigma = 0.001),
    limits  = range(v_pos, na.rm = TRUE),
    oob     = squish,
    breaks  = brks,
    labels  = fmt3,
    na.value = "white",
    name = "Connectivity"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      label.position = "right",
      barwidth  = unit(20, "pt"),
      barheight = unit(80, "pt"),
      ticks = TRUE
    )
  ) +
  labs(x = "TO", y = "FROM") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    panel.grid      = element_blank(),
    axis.title      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y     = element_text(face = "bold")
  ) +
  scale_x_discrete(labels = lab_gc) +
  scale_y_discrete(labels = lab_gc)

print(p_hm)

ggsave(
  filename = file.path(out_dir, "InterCluster_SUMW_heatmap.png"),
  plot     = p_hm,
  width    = 8,
  height   = 8,
  dpi      = 300,
  units    = "in",
  bg = "white"
)

library(ggalluvial)

df_long <- as_tibble(M_SUM0, rownames = "A") |>
  pivot_longer(-A, names_to = "B", values_to = "val") |>
  mutate(
    A = factor(A, levels = clusters),
    B = factor(B, levels = clusters)
  ) |>
  filter(as.character(A) != as.character(B),
         is.finite(val), val > 0)

clusters_lab <- unname(cluster_pretty[clusters])
df_long <- df_long %>%
  mutate(
    A_lab = factor(lab_gc(A), levels = clusters_lab),
    B_lab = factor(lab_gc(B), levels = clusters_lab)
  )

cluster_pal <- c(
  "South-Eastern"    = "#4DBBD5",
  "Central-Eastern"  = "#F39B7F",
  "Central-Western"  = "#3C5488",
  "Northern"         = "#BCBD22",
  "Central-Southern" = "#1B9E77"
)

p_sankey <- ggplot(df_long, aes(y = val, axis1 = A_lab, axis2 = B_lab)) +
  geom_alluvium(aes(fill = A_lab), width = 0.25, alpha = 0.9, knot.pos = 0.35) +
  geom_stratum(width = 0.35, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, vjust = -0.4) +
  scale_x_discrete(limits = c("From", "To")) +
  scale_fill_manual(values = cluster_pal, guide = "none") +
  scale_y_continuous(NULL, breaks = NULL, labels = NULL,
                     expand = expansion(mult = c(0.02, 0.05))) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(face = "bold"),
    panel.grid   = element_blank(),
    plot.margin  = margin(10, 20, 10, 10)
  )

print(p_sankey)

ggsave(
  filename = file.path(out_dir, "Anatidae_noNSbias_Sankey_SUMW.png"),
  plot     = p_sankey,
  width    = 10,
  height   = 6,
  dpi      = 300,
  units    = "in",
  bg = "white"
)

ar_df <- df_long %>%
  left_join(zone_centers %>% select(GeoCluster, x, y), by = c("A" = "GeoCluster")) %>%
  rename(x0 = x, y0 = y) %>%
  left_join(zone_centers %>% select(GeoCluster, x, y), by = c("B" = "GeoCluster")) %>%
  rename(x1 = x, y1 = y) %>%
  mutate(
    From = lab_gc(A),
    To   = lab_gc(B)
  )

thr <- quantile(ar_df$val, 0.75, na.rm = TRUE)
ar_df2 <- ar_df %>% filter(val >= thr)

p_arrows <- ggplot() +
  geom_sf(data = europe, fill = "honeydew", color = "black", linewidth = 0.2) +
  geom_curve(
    data = ar_df2,
    aes(x = x0, y = y0, xend = x1, yend = y1, linewidth = val),
    curvature = 0.25,
    alpha = 0.75,
    arrow = arrow(length = unit(2.2, "mm"), type = "closed")
  ) +
  geom_point(
    data = zone_centers,
    aes(x = x, y = y),
    size = 3, shape = 21, fill = "white", color = "black", stroke = 0.6
  ) +
  geom_text(
    data = zone_centers,
    aes(x = x, y = y, label = Region),
    nudge_y = 1.5,
    size = 4,
    fontface = "bold"
  ) +
  scale_linewidth_continuous(
    name = "SUMW (top 25%)",
    range = c(0.2, 2.0),
    labels = label_number(accuracy = 0.001)
  ) +
  coord_sf(xlim = c(-25, 41), ylim = c(34, 73), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_arrows)

ggsave(
  filename = file.path(out_dir, "InterCluster_SUMW_arrows_top25.png"),
  plot     = p_arrows,
  width    = 10,
  height   = 7,
  dpi      = 300,
  units    = "in",
  bg = "white"
)

message("• Building backbone network (spring+fall pooled) …")

nodes_mig0 <- birds_grid %>%
  filter(season %in% c("spring", "fall")) %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    richness_mig  = n_distinct(genus),
    taxa_list_mig = list(sort(unique(genus))),
    lat           = mean(decimalLatitude,  na.rm = TRUE),
    lon           = mean(decimalLongitude, na.rm = TRUE),
    GeoCluster    = names(which.max(table(GeoCluster))),
    .groups = "drop"
  )

nodes_mig <- merge(grid_05["grid_id"], nodes_mig0, by = "grid_id", all = FALSE)
stopifnot(nrow(nodes_mig) >= 2)

set.seed(42)
ndf <- nodes_mig %>% st_drop_geometry()
coords_mat <- ndf %>% select(lon, lat) %>% as.matrix()
n_nodes    <- nrow(coords_mat)
k_use      <- max(1, min(20, n_nodes - 1))
knn        <- FNN::get.knn(coords_mat, k = k_use)

edges_idx2 <- tibble(
  src = rep(seq_len(n_nodes), each = ncol(knn$nn.index)),
  dst = as.integer(as.vector(knn$nn.index))
) %>%
  filter(src != dst) %>%
  mutate(a = pmin(src, dst), b = pmax(src, dst)) %>%
  distinct(a, b) %>%
  transmute(src = a, dst = b)

eps <- 1e-6
lat_i2 <- ndf$lat[edges_idx2$src]; lon_i2 <- ndf$lon[edges_idx2$src]
lat_j2 <- ndf$lat[edges_idx2$dst]; lon_j2 <- ndf$lon[edges_idx2$dst]

DIST_ij2 <- geodist::geodist(
  data.frame(x = lon_i2, y = lat_i2),
  data.frame(x = lon_j2, y = lat_j2),
  measure = "geodesic",
  paired  = TRUE
) / 1000

taxa <- ndf$taxa_list_mig
jacc_fun <- function(i, j){
  a <- taxa[[i]]; b <- taxa[[j]]
  if (length(a) == 0L || length(b) == 0L) return(0)
  inter <- length(intersect(a, b)); uni <- length(union(a, b))
  if (uni == 0) 0 else inter/uni
}
J_ij2 <- purrr::map2_dbl(edges_idx2$src, edges_idx2$dst, jacc_fun)
W_ij2 <- J_ij2 / pmax(as.numeric(DIST_ij2), eps)

edges_mig <- tibble(
  source_id = ndf$grid_id[edges_idx2$src],
  target_id = ndf$grid_id[edges_idx2$dst],
  W         = as.numeric(W_ij2),
  jaccard   = J_ij2,
  dist_km   = as.numeric(DIST_ij2)
) %>% filter(jaccard > 0, is.finite(W), W > 0)

g_mig <- igraph::graph_from_data_frame(
  d = edges_mig %>%
    mutate(cost = 1 / pmax(W, 1e-9)) %>%
    select(source_id, target_id, cost),
  directed = FALSE,
  vertices = ndf %>% select(grid_id, GeoCluster, richness_mig)
)

bet_mig <- igraph::betweenness(g_mig, directed = FALSE, weights = igraph::E(g_mig)$cost, normalized = TRUE)
bet_vec <- setNames(as.numeric(bet_mig), igraph::V(g_mig)$name)

nodes_mig <- nodes_mig %>%
  mutate(
    betweenness = bet_vec[grid_id],
    backbone    = !is.na(betweenness) & betweenness > 0.01   # tweak if needed
  )

nd_key <- nodes_mig %>%
  st_drop_geometry() %>%
  select(grid_id, GeoCluster, richness_mig, betweenness, backbone)

edges_mig_bb_intra <- edges_mig %>%
  inner_join(nd_key, by = c("source_id" = "grid_id")) %>% rename(Gs = GeoCluster, bb_s = backbone) %>%
  inner_join(nd_key, by = c("target_id" = "grid_id")) %>% rename(Gt = GeoCluster, bb_t = backbone) %>%
  filter(bb_s | bb_t) %>%
  filter(Gs == Gt) %>%
  mutate(GeoCluster = Gs)

nd_xy <- nodes_mig %>% st_drop_geometry() %>% select(grid_id, lon, lat)

edges_bb_intra_sf <- edges_mig_bb_intra %>%
  inner_join(nd_xy, by = c("source_id" = "grid_id")) %>%
  rename(lon_s = lon, lat_s = lat) %>%
  inner_join(nd_xy, by = c("target_id" = "grid_id")) %>%
  rename(lon_t = lon, lat_t = lat) %>%
  mutate(across(c(lon_s, lat_s, lon_t, lat_t), as.numeric)) %>%
  filter(complete.cases(lon_s, lat_s, lon_t, lat_t)) %>%
  mutate(
    geometry = pmap(
      list(lon_s, lat_s, lon_t, lat_t),
      \(xs, ys, xt, yt) st_linestring(matrix(c(xs, ys, xt, yt), ncol = 2, byrow = TRUE))
    )
  ) %>%
  st_as_sf(crs = 4326) %>%
  group_by(GeoCluster) %>%
  mutate(
    W_low  = quantile(W, 0.05, na.rm = TRUE),
    W_high = quantile(W, 0.95, na.rm = TRUE),
    W_clip = pmin(pmax(W, W_low), W_high),
    W_rel  = if_else(W_high > W_low, (W_clip - W_low)/(W_high - W_low), 0.5),
    edge_w = rescale(W_rel, to = c(0.6, 2.2)),
    top    = W >= quantile(W, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()

edges_ctx <- edges_bb_intra_sf
edges_top <- edges_bb_intra_sf %>% filter(top)

nodes_bb_pts <- st_centroid(nodes_mig) %>%
  filter(backbone)

facet_levels <- clusters
nodes_bb_pts$GeoCluster <- factor(nodes_bb_pts$GeoCluster, levels = facet_levels)
edges_ctx$GeoCluster    <- factor(edges_ctx$GeoCluster, levels = facet_levels)
edges_top$GeoCluster    <- factor(edges_top$GeoCluster, levels = facet_levels)

lon_breaks  <- c(-20, 0, 20, 40)
lon_labels  <- c("20°W", "0°", "20°E", "40°E")

p_backbone <- ggplot() +
  geom_sf(data = europe, fill = "honeydew", color = "black", linewidth = 0.2) +
  geom_sf(data = edges_ctx, color = "grey75", linewidth = 0.25, alpha = 0.40) +
  geom_sf(data = edges_top, aes(linewidth = edge_w), color = "brown", alpha = 0.95, lineend = "round") +
  geom_sf(
    data = nodes_bb_pts,
    aes(size = richness_mig, fill = betweenness),
    shape = 21, color = "black", stroke = 0.25, alpha = 0.95
  ) +
  scale_linewidth_continuous(
    name   = "Corridor strength (relative within region)",
    breaks = c(0.7, 1.2, 1.7, 2.2),
    labels = c("weak", "moderate", "strong", "very strong"),
    range  = c(0.6, 2.2)
  ) +
  scale_fill_distiller(
    name      = "Betweenness",
    palette   = "RdYlBu",
    direction = 1,
    limits    = c(0.01, 0.06),
    breaks    = c(0.01, 0.04, 0.06),
    labels    = number_format(accuracy = 0.01)
  ) +
  scale_size_continuous(name = "Migration richness (unique genera)", range = c(1.5, 5)) +
  scale_x_continuous(breaks = lon_breaks, labels = lon_labels) +
  coord_sf(xlim = c(-25, 41), ylim = c(34, 81), expand = FALSE) +
  facet_wrap(
    ~ GeoCluster, ncol = 5,
    labeller = labeller(GeoCluster = as_labeller(cluster_pretty)),
    drop = FALSE
  ) +
  guides(
    fill = guide_colorbar(
      order          = 1,
      title.position = "top",
      direction      = "horizontal",
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.6, "cm")
    ),
    size = guide_legend(
      order = 2,
      title.position = "top"
    ),
    linewidth = guide_legend(
      order = 3,
      title.position = "top",
      override.aes = list(color = "brown")
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "top",
    legend.box         = "horizontal",
    legend.box.spacing = unit(1, "cm"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 11),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(2.0, "lines"),
    panel.spacing.y = unit(2.0, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 2, r = 20, b = 10, l = 10),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(p_backbone)

ggsave(
  filename = file.path(out_dir, "Anatidae_Backbone_noNSbias.png"),
  plot     = p_backbone,
  width    = 14,
  height   = 10,
  dpi      = 300,
  units    = "in",
  bg = "white"
)