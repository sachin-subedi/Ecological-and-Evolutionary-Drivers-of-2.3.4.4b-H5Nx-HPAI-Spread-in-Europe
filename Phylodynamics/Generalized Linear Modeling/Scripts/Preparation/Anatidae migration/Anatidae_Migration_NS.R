suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
  library(FNN)
  library(geodist)
  library(lubridate)
  library(readr)
})

sf_use_s2(TRUE)

csv_path <- "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/output/file/Family/Anatidae/filtered_data_combined_Anatidae.csv"
out_dir  <- "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae_migratory_noradius"
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
k_use   <- max(1, min(20, n_nodes - 1))   # k = 20
knn     <- FNN::get.knn(coords_mat, k = k_use)

edges_idx <- tibble(
  src = rep(seq_len(n_nodes), each = ncol(knn$nn.index)),
  dst = as.integer(as.vector(knn$nn.index))
) %>%
  filter(src != dst) %>%
  mutate(a = pmin(src, dst), b = pmax(src, dst)) %>%
  distinct(a, b) %>%
  transmute(src = a, dst = b)

eps <- 1e-6
lat_i <- nodes_tbl$lat[edges_idx$src]; lon_i <- nodes_tbl$lon[edges_idx$src]
lat_j <- nodes_tbl$lat[edges_idx$dst]; lon_j <- nodes_tbl$lon[edges_idx$dst]
ANGLE_ij <- atan(abs(lat_i - lat_j) / (abs(lon_i - lon_j) + eps))
DIST_ij  <- geodist::geodist(
  data.frame(x = lon_i, y = lat_i),
  data.frame(x = lon_j, y = lat_j),
  measure = "geodesic", paired = TRUE
) / 1000

A_ids <- nodes_tbl$grid_id[edges_idx$src]
B_ids <- nodes_tbl$grid_id[edges_idx$dst]
Gs_static <- nodes_tbl$GeoCluster[edges_idx$src]
Gt_static <- nodes_tbl$GeoCluster[edges_idx$dst]

node_genus_by_year_season <- function(y, s){
  bys <- birds_grid %>% filter(year == y, season == s) %>% st_drop_geometry()
  if (nrow(bys) == 0) return(NULL)
  bys %>%
    group_by(grid_id) %>%
    summarise(
      genus_list = list(sort(unique(genus))),
      abund      = sum(individualCount, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    inner_join(nodes_tbl, by = "grid_id")
}

build_migratory_for_year <- function(y){
  ng_spr <- node_genus_by_year_season(y, "spring")
  spring_df <- NULL
  if (!is.null(ng_spr)) {
    gl <- setNames(ng_spr$genus_list, ng_spr$grid_id)
    cl <- setNames(ng_spr$GeoCluster, ng_spr$grid_id)
    lt <- setNames(ng_spr$lat,        ng_spr$grid_id)
    keep <- (A_ids %in% ng_spr$grid_id) & (B_ids %in% ng_spr$grid_id)
    if (any(keep)) {
      src <- A_ids[keep]; dst <- B_ids[keep]
      lat_i_y <- lt[src];  lat_j_y <- lt[dst]
      Gs <- cl[src];       Gt <- cl[dst]
      jacc <- mapply(function(i,j){
        gi <- gl[[i]]; gj <- gl[[j]]
        if (is.null(gi) || is.null(gj)) 0 else {
          u <- union(gi, gj); if (length(u)==0) 0 else length(intersect(gi,gj))/length(u)
        }
      }, src, dst, USE.NAMES = FALSE)
      keep2 <- jacc > 0
      if (any(keep2)) {
        W <- (as.numeric(ANGLE_ij[keep][keep2]) * jacc[keep2]) / pmax(as.numeric(DIST_ij[keep][keep2]), 1e-6)
        A <- ifelse(lat_i_y[keep2] <= lat_j_y[keep2], Gs[keep2], Gt[keep2])
        B <- ifelse(lat_i_y[keep2] <= lat_j_y[keep2], Gt[keep2], Gs[keep2])
        spring_df <- tibble(A=A, B=B, W=W)
      }
    }
  }
  ng_fall <- node_genus_by_year_season(y, "fall")
  fall_df <- NULL
  if (!is.null(ng_fall)) {
    gl <- setNames(ng_fall$genus_list, ng_fall$grid_id)
    cl <- setNames(ng_fall$GeoCluster, ng_fall$grid_id)
    lt <- setNames(ng_fall$lat,        ng_fall$grid_id)
    keep <- (A_ids %in% ng_fall$grid_id) & (B_ids %in% ng_fall$grid_id)
    if (any(keep)) {
      src <- A_ids[keep]; dst <- B_ids[keep]
      lat_i_y <- lt[src];  lat_j_y <- lt[dst]
      Gs <- cl[src];       Gt <- cl[dst]
      jacc <- mapply(function(i,j){
        gi <- gl[[i]]; gj <- gl[[j]]
        if (is.null(gi) || is.null(gj)) 0 else {
          u <- union(gi, gj); if (length(u)==0) 0 else length(intersect(gi,gj))/length(u)
        }
      }, src, dst, USE.NAMES = FALSE)
      keep2 <- jacc > 0
      if (any(keep2)) {
        W <- (as.numeric(ANGLE_ij[keep][keep2]) * jacc[keep2]) / pmax(as.numeric(DIST_ij[keep][keep2]), 1e-6)
        A <- ifelse(lat_i_y[keep2] >= lat_j_y[keep2], Gs[keep2], Gt[keep2])
        B <- ifelse(lat_i_y[keep2] >= lat_j_y[keep2], Gt[keep2], Gs[keep2])
        fall_df <- tibble(A=A, B=B, W=W)
      }
    }
  }
  mig <- bind_rows(spring_df, fall_df)
  if (is.null(mig) || nrow(mig)==0) return(NULL)
  mig %>%
    filter(A != B) %>%
    group_by(A,B) %>%
    summarise(
      EDGE_COUNT = n(),
      SUMW       = sum(W, na.rm = TRUE),
      MEANW      = mean(W, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = y, season = "migratory")
}

build_nonmig_for_year <- function(y){
  ng <- node_genus_by_year_season(y, "nonmig")
  if (is.null(ng)) return(NULL)
  gl <- setNames(ng$genus_list, ng$grid_id)
  cl <- setNames(ng$GeoCluster, ng$grid_id)
  keep <- (A_ids %in% ng$grid_id) & (B_ids %in% ng$grid_id)
  if (!any(keep)) return(NULL)
  src <- A_ids[keep]; dst <- B_ids[keep]
  Gs <- cl[src];       Gt <- cl[dst]
  jacc <- mapply(function(i,j){
    gi <- gl[[i]]; gj <- gl[[j]]
    if (is.null(gi) || is.null(gj)) 0 else {
      u <- union(gi, gj); if (length(u)==0) 0 else length(intersect(gi,gj))/length(u)
    }
  }, src, dst, USE.NAMES = FALSE)
  keep2 <- jacc > 0
  if (!any(keep2)) return(NULL)
  W <- (as.numeric(ANGLE_ij[keep][keep2]) * jacc[keep2]) / pmax(as.numeric(DIST_ij[keep][keep2]), 1e-6)
  df1 <- tibble(A = Gs[keep2], B = Gt[keep2], W = W)
  df2 <- tibble(A = Gt[keep2], B = Gs[keep2], W = W)
  nm  <- bind_rows(df1, df2) %>%
    filter(A != B) %>%
    group_by(A,B) %>%
    summarise(
      EDGE_COUNT = n(),
      SUMW       = sum(W, na.rm = TRUE),
      MEANW      = mean(W, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = y, season = "nonmig")
  nm
}

message("• Building per-year summaries …")
per_year <- list()
for (y in years) {
  m  <- build_migratory_for_year(y)
  nm <- build_nonmig_for_year(y)
  if (!is.null(m))  per_year[[length(per_year)+1]] <- m
  if (!is.null(nm)) per_year[[length(per_year)+1]] <- nm
}
pair_year <- bind_rows(per_year)
if (nrow(pair_year) == 0) stop("No per-year summaries produced. Check input/date parsing.")
avg_by_season <- function(s, col){
  grid_yab <- tidyr::expand_grid(year = years, A = clusters, B = clusters) %>% filter(A != B)
  dat <- pair_year %>%
    filter(season == s) %>%
    right_join(grid_yab, by = c("year","A","B")) %>%
    mutate(
      EDGE_COUNT = replace_na(EDGE_COUNT, 0),
      SUMW       = replace_na(SUMW, 0),
      MEANW      = replace_na(MEANW, 0)
    )
  dat %>%
    group_by(A,B) %>%
    summarise(val = mean(.data[[col]], na.rm = TRUE), .groups = "drop")
}

ensure_dir(out_dir)

mig_dir <- file.path(out_dir, "season_migratory_raw")
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

message("✓ Wrote MIGRATORY matrices (RAW) to: ", mig_dir)

nm_dir <- file.path(out_dir, "season_nonmig_raw")
ensure_dir(nm_dir)

den_nm <- avg_by_season("nonmig", "EDGE_COUNT")
sum_nm <- avg_by_season("nonmig", "SUMW")
men_nm <- avg_by_season("nonmig", "MEANW")

NM_DEN <- to_square(den_nm, diag_eps = DIAG_EPS)
NM_SUM <- to_square(sum_nm, diag_eps = DIAG_EPS)
NM_MEN <- to_square(men_nm, diag_eps = DIAG_EPS)

write_matrix_tsv(NM_DEN, file.path(nm_dir, "Anatidae_nonmig_EDGE_COUNT.tsv"))
write_matrix_tsv(NM_SUM, file.path(nm_dir, "Anatidae_nonmig_SUMW.tsv"))
write_matrix_tsv(NM_MEN, file.path(nm_dir, "Anatidae_nonmig_MEANW.tsv"))

message("✓ Wrote NONMIG matrices (RAW) to: ", nm_dir)

message("ALL DONE ✅ (Migratory = spring+fall pooled; no radius; RAW only)")


suppressPackageStartupMessages({ library(tidyverse); library(scales); library(grid) })

sum_path_ns <- file.path(out_dir, "season_migratory_raw", "Anatidae_migratory_SUMW.tsv")
M_SUM_NS <- readr::read_tsv(sum_path_ns, show_col_types = FALSE) |>
  column_to_rownames("A") |>
  as.matrix()

stopifnot(all(rownames(M_SUM_NS) %in% clusters), all(colnames(M_SUM_NS) %in% clusters))
M_SUM_NS <- M_SUM_NS[clusters, clusters, drop = FALSE]

hm_df <- as_tibble(M_SUM_NS, rownames = "A") |>
  pivot_longer(-A, names_to = "B", values_to = "val") |>
  mutate(
    A = factor(A, levels = clusters),
    B = factor(B, levels = clusters),
    val_plot = if_else(as.character(A) == as.character(B), NA_real_, val)
  )

v_pos <- hm_df$val_plot[is.finite(hm_df$val_plot) & hm_df$val_plot > 0]
low  <- min(v_pos, na.rm = TRUE)
high <- max(v_pos, na.rm = TRUE)
mid  <- exp(mean(log(v_pos)))
brks <- c(low, mid, high)

fmt3 <- function(x){
  top <- max(v_pos, na.rm = TRUE)
  acc <- if (top < 0.1) 0.001 else if (top < 1) 0.01 else 0.1
  label_number(accuracy = acc)(x)
}

pal_low  <- "#e0f3f8"
pal_mid  <- "#abd9e9"
pal_high <- "#2166ac"

p_hm <- ggplot(hm_df, aes(B, A, fill = val_plot)) +
  geom_tile(color = "black", linewidth = 0.6) +
  scale_fill_gradientn(
    colours = c(pal_low, pal_mid, pal_high),
    trans   = pseudo_log_trans(sigma = 0.001),
    limits  = range(v_pos, na.rm = TRUE),
    oob     = squish,
    breaks  = brks,
    labels  = fmt3,
    na.value = "white",
    name    = "Connectivity"
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
  labs(
    x = "TO",
    y = "FROM",
  ) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    panel.grid      = element_blank(),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y     = element_text(angle = 0, hjust = 1, face = "bold")
  ) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x))

print(p_hm)


ggsave(
  filename = "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae_migratory_noradius/InterCluster_With_NS.png",
  plot     = p_hm,
  width    = 8,
  height   = 8,
  dpi      = 300,
  units    = "in",
  bg = 'white'
)

##Sankey

suppressPackageStartupMessages({ library(tidyverse); library(scales); library(grid) })

sum_path_ns <- file.path(out_dir, "season_migratory_raw", "Anatidae_migratory_SUMW.tsv")
M_SUM_NS <- readr::read_tsv(sum_path_ns, show_col_types = FALSE) |>
  column_to_rownames("A") |>
  as.matrix()

stopifnot(all(rownames(M_SUM_NS) %in% clusters), all(colnames(M_SUM_NS) %in% clusters))
M_SUM_NS <- M_SUM_NS[clusters, clusters, drop = FALSE]

df_long <- as_tibble(M_SUM_NS, rownames = "A") |>
  pivot_longer(-A, names_to = "B", values_to = "val") |>
  mutate(
    A = factor(A, levels = clusters),
    B = factor(B, levels = clusters)
  ) |>
  filter(as.character(A) != as.character(B),
         is.finite(val), val > 0)

clusters_lab <- gsub("_", " ", clusters)
df_long <- df_long %>%
  mutate(
    A_lab = factor(gsub("_", " ", as.character(A)), levels = clusters_lab),
    B_lab = factor(gsub("_", " ", as.character(B)), levels = clusters_lab)
  )

cluster_pal <- c(
  "GeoCluster One"   = "#4DBBD5",  # GeoCluster_One
  "GeoCluster Two"   = "#F39B7F",  # GeoCluster_Two
  "GeoCluster Three" = "#3C5488",  # GeoCluster_Three
  "GeoCluster Four"  = "#BCBD22",  # GeoCluster_Four
  "GeoCluster Five"  = "#1B9E77"   # GeoCluster_Five
)

p_sankey <- ggplot(
  df_long,
  aes(y = val, axis1 = A_lab, axis2 = B_lab)
) +
  geom_alluvium(aes(fill = A_lab), width = 0.25, alpha = 0.9, knot.pos = 0.35) +
  geom_stratum(width = 0.35, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, vjust = -0.4) +
  scale_x_discrete(limits = c("From", "To")) +
  scale_fill_manual(values = cluster_pal, guide = "none") +  # <- your colors
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
  filename = "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae_migratory_noradius/Anatidae_NS_bias_Sankey.png",
  plot     = p_sankey,
  width    = 10,
  height   = 6,
  dpi      = 300,
  units    = "in",
  bg = 'white'
)
###Plotting:

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(FNN); library(geodist); library(igraph); library(purrr)
})

nodes_mig0 <- birds_grid %>%
  dplyr::filter(season %in% c("spring","fall")) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(grid_id) %>%
  dplyr::summarise(
    richness_mig = dplyr::n_distinct(genus),
    taxa_list_mig = list(sort(unique(genus))),
    lat = mean(decimalLatitude,  na.rm = TRUE),
    lon = mean(decimalLongitude, na.rm = TRUE),
    GeoCluster = dplyr::first(GeoCluster),
    .groups = "drop"
  )

if (exists("grid_05")) {
  nodes_mig <- merge(grid_05["grid_id"], nodes_mig0, by = "grid_id", all = FALSE)
} else {
  nodes_mig <- nodes_mig0
}

stopifnot(nrow(nodes_mig) >= 2)

set.seed(42)
coords_mat <- nodes_mig %>% sf::st_drop_geometry() %>% dplyr::select(lon, lat) %>% as.matrix()
n_nodes    <- nrow(coords_mat)
k_use      <- max(1, min(20, n_nodes - 1))
knn        <- FNN::get.knn(coords_mat, k = k_use)

edges_idx <- tibble::tibble(
  src = rep(seq_len(n_nodes), each = ncol(knn$nn.index)),
  dst = as.integer(as.vector(knn$nn.index))
) %>%
  dplyr::filter(src != dst) %>%
  dplyr::mutate(a = pmin(src, dst), b = pmax(src, dst)) %>%
  dplyr::distinct(a, b) %>%
  dplyr::transmute(src = a, dst = b)

eps <- 1e-6
ndf <- nodes_mig %>% sf::st_drop_geometry()

lat_i <- ndf$lat[edges_idx$src]; lon_i <- ndf$lon[edges_idx$src]
lat_j <- ndf$lat[edges_idx$dst]; lon_j <- ndf$lon[edges_idx$dst]

ANGLE_ij <- atan(abs(lat_i - lat_j) / (abs(lon_i - lon_j) + eps))
DIST_ij  <- geodist::geodist(
  data.frame(x = lon_i, y = lat_i),
  data.frame(x = lon_j, y = lat_j),
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
J_ij <- purrr::map2_dbl(edges_idx$src, edges_idx$dst, jacc_fun)

W_ij <- (as.numeric(ANGLE_ij) * J_ij) / pmax(as.numeric(DIST_ij), eps)
edges_mig <- tibble::tibble(
  source_id = ndf$grid_id[edges_idx$src],
  target_id = ndf$grid_id[edges_idx$dst],
  W         = W_ij,
  jaccard   = J_ij,
  dist_km   = as.numeric(DIST_ij),
  angle     = as.numeric(ANGLE_ij)
) %>% dplyr::filter(jaccard > 0, is.finite(W), W > 0)

g_mig <- igraph::graph_from_data_frame(
  d = edges_mig %>% dplyr::mutate(cost = 1 / pmax(W, 1e-9)) %>% dplyr::select(source_id, target_id, cost),
  directed = FALSE,
  vertices = ndf %>% dplyr::select(grid_id, GeoCluster, richness_mig)
)
bet_mig <- igraph::betweenness(g_mig, directed = FALSE, weights = igraph::E(g_mig)$cost, normalized = TRUE)
bet_vec <- setNames(as.numeric(bet_mig), igraph::V(g_mig)$name)

if ("geometry" %in% names(nodes_mig)) {
  nodes_mig <- nodes_mig %>%
    dplyr::mutate(
      betweenness = bet_vec[grid_id],
      backbone    = !is.na(betweenness) & betweenness > 0.01   # or cluster-wise top quartile
    )
} else {
  nodes_mig <- nodes_mig %>%
    dplyr::mutate(
      betweenness = bet_vec[grid_id],
      backbone    = !is.na(betweenness) & betweenness > 0.01
    )
}

nd_key <- nodes_mig %>% sf::st_drop_geometry() %>% dplyr::select(grid_id, GeoCluster, richness_mig, betweenness, backbone)
edges_mig_bb_intra <- edges_mig %>%
  dplyr::inner_join(nd_key, by = c("source_id" = "grid_id")) %>% dplyr::rename(Gs = GeoCluster) %>%
  dplyr::inner_join(nd_key, by = c("target_id" = "grid_id")) %>% dplyr::rename(Gt = GeoCluster) %>%
  dplyr::filter(backbone.x | backbone.y) %>%
  dplyr::filter(Gs == Gt) %>%
  dplyr::mutate(GeoCluster = Gs)



library(dplyr)
library(sf)
library(purrr)
library(ggplot2)
library(scales)

if (!("geometry" %in% names(nodes_mig))) {
  nodes_mig_pts <- st_as_sf(nodes_mig, coords = c("lon","lat"), crs = 4326)
} else {
  nodes_mig_pts <- st_centroid(nodes_mig)
}

nodes_bb <- nodes_mig %>% mutate(has_geom = "geometry" %in% names(.)) %>%
  { if (.$has_geom[1]) . else st_as_sf(., coords = c("lon","lat"), crs = 4326) } %>%
  filter(backbone)

nodes_bb_pts <- if ("geometry" %in% names(nodes_bb)) st_centroid(nodes_bb) else nodes_bb
nodes_bb_pts$GeoCluster      <- factor(nodes_bb_pts$GeoCluster,
                                       levels = c("GeoCluster_One","GeoCluster_Two","GeoCluster_Three","GeoCluster_Four","GeoCluster_Five")
)

nd_xy <- nodes_mig %>%
  st_drop_geometry() %>%
  select(grid_id, lon, lat)

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
  mutate(edge_w = rescale(W, to = c(0.4, 2.0)))

# --- 3) Plot -------------------------------------------------------------------
ggplot() +
  geom_sf(data = europe, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = edges_bb_intra_sf, aes(linewidth = edge_w), color = "grey25", alpha = 0.85) +
  geom_sf(
    data = nodes_bb_pts,
    aes(size = richness_mig, fill = betweenness),
    shape = 21, color = "black", alpha = 0.95
  ) +
  scale_linewidth_continuous(name = "Corridor strength (W)") +
  scale_fill_viridis_c(name = "Betweenness") +
  scale_size_continuous(name = "Migration richness (unique genera)", range = c(1.4, 6)) +
  coord_sf(xlim = c(-15, 40), ylim = c(35, 80), expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right") +
  labs(
    title = "Anatidae — Migration Backbone (within cluster)",
    subtitle = "Lines ∝ W; node size = migration-season genus richness"
  ) +
  facet_wrap(~ GeoCluster)

library(dplyr)
library(sf)
library(ggplot2)
library(scales)

facet_levels <- c("GeoCluster_One","GeoCluster_Two","GeoCluster_Three",
                  "GeoCluster_Four","GeoCluster_Five")

nodes_bb_pts$GeoCluster      <- factor(nodes_bb_pts$GeoCluster,      levels = facet_levels)
edges_ctx$GeoCluster         <- factor(edges_ctx$GeoCluster,         levels = facet_levels)
edges_top$GeoCluster         <- factor(edges_top$GeoCluster,         levels = facet_levels)

nice_names <- c(
  GeoCluster_One   = "GeoCluster One",
  GeoCluster_Two   = "GeoCluster Two",
  GeoCluster_Three = "GeoCluster Three",
  GeoCluster_Four  = "GeoCluster Four",
  GeoCluster_Five  = "GeoCluster Five"
)

edges_bb_intra_sf <- edges_bb_intra_sf %>%
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

edges_ctx  <- dplyr::filter(edges_bb_intra_sf,  TRUE)
edges_top  <- dplyr::filter(edges_bb_intra_sf, top)



library(grid)
lon_breaks  <- c(-20, 0, 20, 40)
lon_labels  <- c("20°W", "0°", "20°E", "40°E")

p <- ggplot() +
  geom_sf(data = europe, fill = "honeydew", color = "black", linewidth = 0.2) +
  geom_sf(data = edges_ctx, color = "grey75", linewidth = 0.25, alpha = 0.40) +
  geom_sf(
    data = edges_top,
    aes(linewidth = edge_w),
    color = "brown",
    alpha = 0.95,
    lineend = "round"
  ) +
  geom_sf(
    data = nodes_bb_pts,
    aes(size = richness_mig, fill = betweenness),
    shape = 21, color = "black", stroke = 0.25, alpha = 0.95
  ) +
  scale_linewidth_continuous(
    name   = "Corridor strength (relative within cluster)",
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
  scale_size_continuous(
    name  = "Migration richness (unique genera)",
    range = c(1.5, 5)
  ) +
  scale_x_continuous(breaks = lon_breaks, labels = lon_labels) +
  coord_sf(xlim = c(-25, 41), ylim = c(34, 81), expand = FALSE) +
  facet_wrap(
    ~ GeoCluster, ncol = 5,
    labeller = labeller(GeoCluster = as_labeller(nice_names)),
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
    text = element_text(size = 14),
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 11),

    legend.position    = "top",
    legend.box         = "horizontal",
    legend.box.spacing = unit(1, "cm"),
    
    panel.grid.major = element_line(color = "grey92", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    
    panel.spacing.x = unit(2.0, "lines"),
    panel.spacing.y = unit(2.0, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = margin(t = 2, r = 20, b = 10, l = 10)
  )
p

ggsave(
  filename = "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae_migratory_noradius/Anatidae_NS.png",
  plot     = p,
  width    = 14,
  height   = 10,
  dpi      = 300,
  units    = "in",
  bg = 'white'
)




library(dplyr)
library(readr)
mig <- birds_grid %>%
  dplyr::filter(season %in% c("spring","fall")) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(grid_id, GeoCluster, genus)

n_cells_total   <- mig %>% distinct(grid_id) %>% nrow()
n_clusters_total<- mig %>% distinct(GeoCluster) %>% nrow()


common_overall <- mig %>%
  distinct(grid_id, genus) %>%
  count(genus, name = "n_cells") %>%
  mutate(prop_cells = n_cells / n_cells_total) %>%
  arrange(desc(n_cells))

top25_overall <- common_overall %>% slice_max(n_cells, n = 25)
print(top25_overall, n = 25)

write_tsv(common_overall,
          "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae/common_genera_overall.tsv")


library(dplyr)

mig_ccg <- mig %>%
  distinct(GeoCluster, grid_id, genus)

cluster_cells <- mig_ccg %>%
  distinct(GeoCluster, grid_id) %>%
  count(GeoCluster, name = "cells_in_cluster")

common_by_cluster <- mig_ccg %>%
  count(GeoCluster, genus, name = "n_cells") %>%
  left_join(cluster_cells, by = "GeoCluster") %>%
  mutate(
    prop_cells = n_cells / cells_in_cluster
  ) %>%
  arrange(GeoCluster, desc(n_cells))

top10_by_cluster <- common_by_cluster %>%
  group_by(GeoCluster) %>%
  slice_max(n_cells, n = 10, with_ties = FALSE) %>%
  ungroup()
top10_by_cluster 
readr::write_tsv(common_by_cluster, "/Users/sachinsubedi/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/eBird/predictors/Anatidae/common_genera_by_geocluster.tsv")

coverage_thresh <- 0.40
core_coverage <- common_overall %>%
  filter(prop_cells >= coverage_thresh) %>%
  arrange(desc(prop_cells))

pancluster_core <- mig %>%
  distinct(genus, GeoCluster) %>%
  count(genus, name = "n_clusters") %>%
  mutate(prop_clusters = n_clusters / n_clusters_total) %>%
  arrange(desc(n_clusters))

core_all5   <- pancluster_core %>% filter(n_clusters == n_clusters_total)
core_4of5   <- pancluster_core %>% filter(n_clusters >= n_clusters_total - 1)

print(core_coverage, n = 50)
print(core_all5,   n = 50)


