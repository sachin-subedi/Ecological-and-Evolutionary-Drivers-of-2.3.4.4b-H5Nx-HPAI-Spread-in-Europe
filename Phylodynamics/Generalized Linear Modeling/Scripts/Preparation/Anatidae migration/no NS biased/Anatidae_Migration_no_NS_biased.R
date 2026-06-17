suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
  library(FNN),
  library(geodist)
  library(lubridate)
  library(readr)
  library(ggalluvial)
  library(igraph)
  library(purrr)
  library(scales)
})

sf_use_s2(TRUE)

csv_path <- "filtered_data_combined_Anatidae.csv"
out_dir  <- "no_NS_biased/"

years    <- 2016:2025
DIAG_EPS <- 1e-11

ensure_dir <- function(p){
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
ensure_dir(out_dir)

jacc_fast <- function(li, lj){
  if (is.null(li) || is.null(lj)) return(0)
  inter <- length(intersect(li, lj)); uni <- length(union(li, lj))
  if (uni == 0) 0 else inter / uni
}

season_from_doy <- function(d){
  d <- as.integer(d)
  ifelse(d >= 46  & d < 135, "spring",
         ifelse(d >= 258 & d < 319, "fall", "nonmig"))
}

clusters_10 <- c("HC_Cluster_1_Alpine",
                 "HC_Cluster_1_Atlantic",
                 "HC_Cluster_1_Continental",
                 "HC_Cluster_2_Alpine",
                 "HC_Cluster_2_Continental",
                 "HC_Cluster_2_Mediterranean",
                 "HC_Cluster_2_Pannonian",
                 "HC_Cluster_3_Alpine",
                 "HC_Cluster_3_Boreal",
                 "HC_Cluster_4_Mediterranean")

clusters <- clusters_10

cluster_pretty_10 <- c(
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

cluster_pal_10 <- c(
  "HC_Cluster_1_Alpine"        = "#0072B2",
  "HC_Cluster_1_Atlantic"      = "#CD5C5C",
  "HC_Cluster_1_Continental"   = "#BCBD22",
  "HC_Cluster_2_Alpine"        = "#26A69A",
  "HC_Cluster_2_Continental"   = "#CC79A7",
  "HC_Cluster_2_Mediterranean" = "lightslategray",
  "HC_Cluster_2_Pannonian"     = "#6B4C1B",
  "HC_Cluster_3_Alpine"        = "#7570B3",
  "HC_Cluster_3_Boreal"        = "#006D2C",
  "HC_Cluster_4_Mediterranean" = "#56B4E9"
)

lab_gc10 <- function(x){
  x <- as.character(x)
  out <- unname(cluster_pretty_10[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

clusters_4 <- c("HC_Cluster_1_Atlantic",
                "HC_Cluster_1_Continental",
                "HC_Cluster_2_Continental",
                "HC_Cluster_3_Boreal")

cluster_pretty_4 <- c(
  "HC_Cluster_1_Atlantic"    = "Atlantic",
  "HC_Cluster_1_Continental" = "Western Continental",
  "HC_Cluster_2_Continental" = "Eastern Continental",
  "HC_Cluster_3_Boreal"      = "Boreal Baltic"
)

cluster_pal_4 <- c(
  "HC_Cluster_1_Atlantic"    = "#CD5C5C",
  "HC_Cluster_1_Continental" = "#BCBD22",
  "HC_Cluster_2_Continental" = "#0072B2",
  "HC_Cluster_3_Boreal"      = "#006D2C"
)

lab_gc4 <- function(x){
  x <- as.character(x)
  out <- unname(cluster_pretty_4[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

cluster_region_lookup <- tribble(
  ~Country,                    ~GeoCluster,
  "Albania",                   "HC_Cluster_2_Mediterranean",
  "Austria",                   "HC_Cluster_2_Alpine",
  "Belgium",                   "HC_Cluster_1_Atlantic",
  "Bosnia and Herzegovina",    "HC_Cluster_2_Alpine",
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
  "Netherlands",               "HC_Cluster_1_Atlantic",
  "North Macedonia",           "HC_Cluster_2_Continental",
  "Norway",                    "HC_Cluster_3_Alpine",
  "Poland",                    "HC_Cluster_2_Continental",
  "Portugal",                  "HC_Cluster_4_Mediterranean",
  "Romania",                   "HC_Cluster_2_Continental",
  "Republic of Serbia",        "HC_Cluster_2_Continental",
  "Slovakia",                  "HC_Cluster_2_Alpine",
  "Slovenia",                  "HC_Cluster_2_Continental",
  "Spain",                     "HC_Cluster_4_Mediterranean",
  "Sweden",                    "HC_Cluster_3_Boreal",
  "Switzerland",               "HC_Cluster_1_Alpine",
  "Ukraine",                   "HC_Cluster_2_Continental",
  "United Kingdom",            "HC_Cluster_1_Atlantic"
) %>% as.data.frame()

to_square <- function(df_AB_val, diag_eps = DIAG_EPS){
  grid <- tidyr::expand_grid(A = clusters, B = clusters)
  mdf  <- grid %>%
    left_join(df_AB_val, by = c("A","B")) %>%
    mutate(val = replace_na(val, 0),
           A   = factor(A, levels = clusters),
           B   = factor(B, levels = clusters)) %>%
    arrange(A, B)
  wide <- mdf %>% select(A, B, val) %>%
    pivot_wider(names_from = B, values_from = val) %>% arrange(A)
  mm <- as.matrix(wide[, -1, drop = FALSE]); rownames(mm) <- wide$A
  diag(mm) <- diag_eps
  mm
}

write_matrix_tsv <- function(mat, file_path){
  out <- tibble(A = rownames(mat)) %>%
    bind_cols(as_tibble(mat, .name_repair = "minimal"))
  readr::write_tsv(out, file_path)
}

message("• Loading bird data …")
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

europe_full <- ne_countries(continent = "Europe", returnclass = "sf") %>%
  filter(name != "Russia") %>%
  st_make_valid()

europe <- europe_full %>%
  transmute(Country = name) %>%
  st_make_valid()

europe_gc <- europe
europe_gc$GeoCluster <- cluster_region_lookup$GeoCluster[
  match(europe_gc$Country, cluster_region_lookup$Country)
]
europe_gc <- europe_gc[!is.na(europe_gc$GeoCluster), ]

birds <- st_join(birds, europe_gc, join = st_within, left = FALSE)
birds <- birds[!is.na(birds$GeoCluster), ]

cols_b <- names(birds)
if (all(c("year","month","day") %in% cols_b)){
  birds$year  <- as.integer(birds$year)
  birds$month <- pmax(1L, pmin(12L, as.integer(birds$month)))
  d_raw       <- suppressWarnings(as.integer(birds$day))
  birds$day   <- ifelse(is.na(d_raw) | d_raw < 1L | d_raw > 31L, 1L, d_raw)
  birds$date  <- make_date(year = birds$year, month = birds$month, day = birds$day)
} else if (all(c("year","startDayOfYear") %in% cols_b)){
  birds$year           <- as.integer(birds$year)
  birds$startDayOfYear <- pmax(1, pmin(366, as.integer(round(birds$startDayOfYear))))
  birds$date           <- ymd(paste0(birds$year,"-01-01")) + days(birds$startDayOfYear - 1)
} else stop("No usable date fields.")

birds$doy    <- yday(birds$date)
birds$season <- season_from_doy(birds$doy)
birds        <- birds %>% filter(year %in% years)

eu_3035       <- st_transform(europe, 3035)
eu_union_3035 <- st_union(eu_3035)
europe_union  <- st_transform(eu_union_3035, 4326)

grid_05 <- st_make_grid(europe_union, cellsize = c(0.5,0.5),
                        what = "polygons", square = TRUE) %>%
  st_sf(crs = 4326) %>%
  st_filter(europe_union, .predicate = st_intersects)

cc <- st_coordinates(st_centroid(st_geometry(grid_05)))
grid_05$grid_id <- paste0("N_", round(cc[,"Y"],3), "_", round(cc[,"X"],3))

birds_grid <- st_join(birds, grid_05["grid_id"], join = st_within, left = FALSE)

nodes_tbl <- birds_grid %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    lat        = mean(decimalLatitude,  na.rm = TRUE),
    lon        = mean(decimalLongitude, na.rm = TRUE),
    Country    = names(which.max(table(Country))),
    GeoCluster = names(which.max(table(GeoCluster))),
    .groups = "drop"
  )

set.seed(42)
coords_all <- nodes_tbl %>% select(lon, lat) %>% as.matrix()
n_all      <- nrow(coords_all)
k_all      <- max(1, min(20, n_all - 1))
knn_all    <- FNN::get.knn(coords_all, k = k_all)

edges_all <- tibble(
  src = rep(seq_len(n_all), each = ncol(knn_all$nn.index)),
  dst = as.integer(as.vector(knn_all$nn.index))
) %>%
  filter(src != dst) %>%
  mutate(a = pmin(src,dst), b = pmax(src,dst)) %>%
  distinct(a, b) %>%
  transmute(src = a, dst = b)

lat_ai <- nodes_tbl$lat[edges_all$src]; lon_ai <- nodes_tbl$lon[edges_all$src]
lat_aj <- nodes_tbl$lat[edges_all$dst]; lon_aj <- nodes_tbl$lon[edges_all$dst]

DIST_all <- geodist::geodist(
  data.frame(x = lon_ai, y = lat_ai),
  data.frame(x = lon_aj, y = lat_aj),
  measure = "geodesic", paired = TRUE
) / 1000

A_ids_all <- nodes_tbl$grid_id[edges_all$src]
B_ids_all <- nodes_tbl$grid_id[edges_all$dst]

node_genus_by_year_season <- function(y, s){
  bys <- birds_grid %>% filter(year == y, season == s) %>% st_drop_geometry()
  if (nrow(bys) == 0) return(NULL)
  bys %>%
    group_by(grid_id) %>%
    summarise(
      genus_list = list(sort(unique(genus))),
      .groups = "drop"
    ) %>%
    inner_join(nodes_tbl, by = "grid_id")
}

build_migratory_for_year <- function(y){

  ng_spr    <- node_genus_by_year_season(y, "spring")
  spring_df <- NULL
  if (!is.null(ng_spr)){
    gl   <- setNames(ng_spr$genus_list, ng_spr$grid_id)
    cl   <- setNames(ng_spr$GeoCluster, ng_spr$grid_id)
    lt   <- setNames(ng_spr$lat,        ng_spr$grid_id)
    keep <- (A_ids_all %in% ng_spr$grid_id) & (B_ids_all %in% ng_spr$grid_id)
    if (any(keep)){
      src <- A_ids_all[keep]; dst <- B_ids_all[keep]
      lat_i_y <- lt[src]; lat_j_y <- lt[dst]
      Gs <- cl[src];      Gt <- cl[dst]
      jacc <- mapply(function(i,j){
        gi <- gl[[i]]; gj <- gl[[j]]
        if (is.null(gi)||is.null(gj)) return(0)
        u <- union(gi,gj); if(length(u)==0) 0 else length(intersect(gi,gj))/length(u)
      }, src, dst, USE.NAMES = FALSE)
      keep2 <- jacc > 0
      if (any(keep2)){
        W <- jacc[keep2] / pmax(as.numeric(DIST_all[keep][keep2]), 1e-6)
        A <- ifelse(lat_i_y[keep2] <= lat_j_y[keep2], Gs[keep2], Gt[keep2])
        B <- ifelse(lat_i_y[keep2] <= lat_j_y[keep2], Gt[keep2], Gs[keep2])
        spring_df <- tibble(A=A, B=B, W=W)
      }
    }
  }
  
  ng_fall  <- node_genus_by_year_season(y, "fall")
  fall_df  <- NULL
  if (!is.null(ng_fall)){
    gl   <- setNames(ng_fall$genus_list, ng_fall$grid_id)
    cl   <- setNames(ng_fall$GeoCluster, ng_fall$grid_id)
    lt   <- setNames(ng_fall$lat,        ng_fall$grid_id)
    keep <- (A_ids_all %in% ng_fall$grid_id) & (B_ids_all %in% ng_fall$grid_id)
    if (any(keep)){
      src <- A_ids_all[keep]; dst <- B_ids_all[keep]
      lat_i_y <- lt[src]; lat_j_y <- lt[dst]
      Gs <- cl[src];      Gt <- cl[dst]
      jacc <- mapply(function(i,j){
        gi <- gl[[i]]; gj <- gl[[j]]
        if (is.null(gi)||is.null(gj)) return(0)
        u <- union(gi,gj); if(length(u)==0) 0 else length(intersect(gi,gj))/length(u)
      }, src, dst, USE.NAMES = FALSE)
      keep2 <- jacc > 0
      if (any(keep2)){
        W <- jacc[keep2] / pmax(as.numeric(DIST_all[keep][keep2]), 1e-6)
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
    summarise(EDGE_COUNT=n(), SUMW=sum(W,na.rm=TRUE), MEANW=mean(W,na.rm=TRUE), .groups="drop") %>%
    mutate(year=y, season="migratory")
}

build_nonmig_for_year <- function(y){
  ng <- node_genus_by_year_season(y, "nonmig")
  if (is.null(ng)) return(NULL)
  gl   <- setNames(ng$genus_list, ng$grid_id)
  cl   <- setNames(ng$GeoCluster, ng$grid_id)
  keep <- (A_ids_all %in% ng$grid_id) & (B_ids_all %in% ng$grid_id)
  if (!any(keep)) return(NULL)
  src <- A_ids_all[keep]; dst <- B_ids_all[keep]
  Gs <- cl[src]; Gt <- cl[dst]
  jacc <- mapply(function(i,j){
    gi <- gl[[i]]; gj <- gl[[j]]
    if (is.null(gi)||is.null(gj)) return(0)
    u <- union(gi,gj); if(length(u)==0) 0 else length(intersect(gi,gj))/length(u)
  }, src, dst, USE.NAMES = FALSE)
  keep2 <- jacc > 0
  if (!any(keep2)) return(NULL)
  W <- jacc[keep2] / pmax(as.numeric(DIST_all[keep][keep2]), 1e-6)
  bind_rows(
    tibble(A=Gs[keep2], B=Gt[keep2], W=W),
    tibble(A=Gt[keep2], B=Gs[keep2], W=W)
  ) %>%
    filter(A != B) %>%
    group_by(A,B) %>%
    summarise(EDGE_COUNT=n(), SUMW=sum(W,na.rm=TRUE), MEANW=mean(W,na.rm=TRUE), .groups="drop") %>%
    mutate(year=y, season="nonmig")
}

message("• Building per-year matrices …")
per_year <- list()
for (y in years){
  m  <- build_migratory_for_year(y)
  nm <- build_nonmig_for_year(y)
  if (!is.null(m))  per_year[[length(per_year)+1]] <- m
  if (!is.null(nm)) per_year[[length(per_year)+1]] <- nm
}
pair_year <- bind_rows(per_year)
if (nrow(pair_year)==0) stop("No per-year summaries. Check input/date parsing.")

avg_by_season <- function(s, col){
  grid_yab <- tidyr::expand_grid(year=years, A=clusters, B=clusters) %>% filter(A!=B)
  pair_year %>%
    filter(season==s) %>%
    right_join(grid_yab, by=c("year","A","B")) %>%
    mutate(EDGE_COUNT=replace_na(EDGE_COUNT,0),
           SUMW=replace_na(SUMW,0),
           MEANW=replace_na(MEANW,0)) %>%
    group_by(A,B) %>%
    summarise(val=mean(.data[[col]], na.rm=TRUE), .groups="drop")
}

mig_dir <- file.path(out_dir, "season_migratory_raw")
ensure_dir(mig_dir)

M_DEN <- to_square(avg_by_season("migratory","EDGE_COUNT"))
M_SUM <- to_square(avg_by_season("migratory","SUMW"))
M_MEN <- to_square(avg_by_season("migratory","MEANW"))

write_matrix_tsv(M_DEN, file.path(mig_dir, "Anatidae_migratory_EDGE_COUNT_No_NS.tsv"))
write_matrix_tsv(M_SUM, file.path(mig_dir, "Anatidae_migratory_SUMW_No_NS.tsv"))
write_matrix_tsv(M_MEN, file.path(mig_dir, "Anatidae_migratory_MEANW_No_NS.tsv"))
message("✓ Wrote MIGRATORY matrices → ", mig_dir)

nm_dir <- file.path(out_dir, "season_nonmig_raw")
ensure_dir(nm_dir)

NM_DEN <- to_square(avg_by_season("nonmig","EDGE_COUNT"))
NM_SUM <- to_square(avg_by_season("nonmig","SUMW"))
NM_MEN <- to_square(avg_by_season("nonmig","MEANW"))

write_matrix_tsv(NM_DEN, file.path(nm_dir, "Anatidae_nonmig_EDGE_COUNT_No_NS.tsv"))
write_matrix_tsv(NM_SUM, file.path(nm_dir, "Anatidae_nonmig_SUMW_No_NS.tsv"))
write_matrix_tsv(NM_MEN, file.path(nm_dir, "Anatidae_nonmig_MEANW_No_NS.tsv"))
message("✓ Wrote NONMIG matrices → ", nm_dir)

M_SUM_NS <- M_SUM

nodes_mig0 <- birds_grid %>%
  filter(season %in% c("spring","fall")) %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    richness_mig  = n_distinct(genus),
    taxa_list_mig = list(sort(unique(genus))),
    lat           = mean(decimalLatitude,  na.rm = TRUE),
    lon           = mean(decimalLongitude, na.rm = TRUE),
    GeoCluster    = names(which.max(table(GeoCluster))),
    .groups       = "drop"
  )

nodes_mig <- merge(grid_05["grid_id"], nodes_mig0, by = "grid_id", all = FALSE)
ndf       <- nodes_mig %>% st_drop_geometry()

# ── kNN edges on migration nodes ─────────────────────────────────────────────
set.seed(42)
coords_m <- ndf %>% select(lon, lat) %>% as.matrix()
n_m      <- nrow(coords_m)
k_m      <- max(1, min(20, n_m - 1))
knn_m    <- FNN::get.knn(coords_m, k = k_m)

edges_idx2 <- tibble(
  src = rep(seq_len(n_m), each = ncol(knn_m$nn.index)),
  dst = as.integer(as.vector(knn_m$nn.index))
) %>%
  filter(src != dst) %>%
  mutate(a = pmin(src,dst), b = pmax(src,dst)) %>%
  distinct(a, b) %>%
  transmute(src = a, dst = b)

DIST_ij2 <- geodist::geodist(
  data.frame(x = ndf$lon[edges_idx2$src], y = ndf$lat[edges_idx2$src]),
  data.frame(x = ndf$lon[edges_idx2$dst], y = ndf$lat[edges_idx2$dst]),
  measure = "geodesic", paired = TRUE
) / 1000

taxa  <- ndf$taxa_list_mig
J_ij2 <- purrr::map2_dbl(edges_idx2$src, edges_idx2$dst,
                         \(i,j) jacc_fast(taxa[[i]], taxa[[j]]))

W_ij2 <- J_ij2 / pmax(as.numeric(DIST_ij2), 1e-6)

edges_mig <- tibble(
  source_id = ndf$grid_id[edges_idx2$src],
  target_id = ndf$grid_id[edges_idx2$dst],
  W         = as.numeric(W_ij2),
  jaccard   = J_ij2,
  dist_km   = as.numeric(DIST_ij2)
) %>% filter(jaccard > 0, is.finite(W), W > 0)

build_backbone <- function(cluster_vec, edges_mig, ndf, nodes_mig, bet_threshold = 0.01){
  
  nodes_sub <- nodes_mig %>% filter(GeoCluster %in% cluster_vec)
  ndf_sub   <- nodes_sub %>% st_drop_geometry()
  
  edges_sub <- edges_mig %>%
    inner_join(ndf_sub %>% select(grid_id, GeoCluster), by = c("source_id"="grid_id")) %>%
    rename(Gs = GeoCluster) %>%
    inner_join(ndf_sub %>% select(grid_id, GeoCluster), by = c("target_id"="grid_id")) %>%
    rename(Gt = GeoCluster) %>%
    filter(Gs %in% cluster_vec, Gt %in% cluster_vec)
  
  if (nrow(edges_sub)==0) return(list(nodes=NULL, edges_ctx=NULL, edges_top=NULL))
  
  g_sub <- igraph::graph_from_data_frame(
    d        = edges_sub %>% mutate(cost=1/pmax(W,1e-9)) %>% select(source_id,target_id,cost),
    directed = FALSE,
    vertices = ndf_sub %>% select(grid_id, GeoCluster, richness_mig)
  )
  bet   <- igraph::betweenness(g_sub, directed=FALSE,
                               weights=igraph::E(g_sub)$cost, normalized=TRUE)
  bet_v <- setNames(as.numeric(bet), igraph::V(g_sub)$name)
  
  nodes_sub <- nodes_sub %>%
    mutate(betweenness = bet_v[grid_id],
           backbone    = !is.na(betweenness) & betweenness > bet_threshold)
  
  nd_key <- nodes_sub %>%
    st_drop_geometry() %>%
    select(grid_id, GeoCluster, richness_mig, betweenness, backbone)
  
  edges_bb <- edges_sub %>%
    inner_join(nd_key, by=c("source_id"="grid_id")) %>% rename(Gs2=GeoCluster, bb_s=backbone) %>%
    inner_join(nd_key, by=c("target_id"="grid_id")) %>% rename(Gt2=GeoCluster, bb_t=backbone) %>%
    filter(bb_s | bb_t, Gs2==Gt2) %>%
    mutate(GeoCluster=Gs2)
  
  nd_xy <- nodes_sub %>% st_drop_geometry() %>% select(grid_id, lon, lat)
  
  edges_sf <- edges_bb %>%
    inner_join(nd_xy, by=c("source_id"="grid_id")) %>% rename(lon_s=lon, lat_s=lat) %>%
    inner_join(nd_xy, by=c("target_id"="grid_id")) %>% rename(lon_t=lon, lat_t=lat) %>%
    mutate(across(c(lon_s,lat_s,lon_t,lat_t), as.numeric)) %>%
    filter(complete.cases(lon_s,lat_s,lon_t,lat_t)) %>%
    mutate(geometry = pmap(list(lon_s,lat_s,lon_t,lat_t),
                           \(xs,ys,xt,yt) st_linestring(matrix(c(xs,ys,xt,yt),ncol=2,byrow=TRUE)))) %>%
    st_as_sf(crs=4326) %>%
    group_by(GeoCluster) %>%
    mutate(
      W_low  = quantile(W, 0.05, na.rm=TRUE),
      W_high = quantile(W, 0.95, na.rm=TRUE),
      W_clip = pmin(pmax(W, W_low), W_high),
      W_rel  = if_else(W_high>W_low, (W_clip-W_low)/(W_high-W_low), 0.5),
      edge_w = rescale(W_rel, to=c(0.6,2.2)),
      top    = W >= quantile(W, 0.75, na.rm=TRUE)
    ) %>%
    ungroup() %>%
    mutate(GeoCluster = factor(GeoCluster, levels=cluster_vec))
  
  nodes_pts <- st_centroid(nodes_sub) %>%
    filter(backbone) %>%
    mutate(GeoCluster = factor(GeoCluster, levels=cluster_vec))
  
  list(nodes=nodes_pts, edges_ctx=edges_sf, edges_top=edges_sf %>% filter(top))
}

message("• Building 10-cluster backbone …")
bb10 <- build_backbone(clusters_10, edges_mig, ndf, nodes_mig)

message("• Building 4-cluster backbone …")
bb4  <- build_backbone(clusters_4,  edges_mig, ndf, nodes_mig)

build_sankey <- function(cluster_vec, pretty_vec, pal_vec, M){
  
  M_sub <- M[cluster_vec, cluster_vec, drop=FALSE]
  labs  <- unname(pretty_vec[cluster_vec])
  
  df <- as_tibble(M_sub, rownames="A") %>%
    pivot_longer(-A, names_to="B", values_to="val") %>%
    filter(A %in% cluster_vec, B %in% cluster_vec) %>%
    mutate(A=factor(A,levels=cluster_vec), B=factor(B,levels=cluster_vec)) %>%
    filter(as.character(A)!=as.character(B), is.finite(val), val>0) %>%
    mutate(
      A_lab  = factor(unname(pretty_vec[as.character(A)]), levels=labs),
      B_lab  = factor(unname(pretty_vec[as.character(B)]), levels=labs),
      A_lab2 = paste0(A_lab,"_from"),
      B_lab2 = paste0(B_lab,"_to")
    )
  
  fill_vals <- setNames(unname(pal_vec), unname(pretty_vec[names(pal_vec)]))
  
  ggplot(df, aes(y=val, axis1=A_lab2, axis2=B_lab2)) +
    geom_alluvium(aes(fill=A_lab), width=0.2, alpha=0.85, knot.pos=0.4) +
    geom_stratum(width=0.3, fill="white", color="black") +
    geom_text(stat="stratum",
              aes(label=gsub("_(from|to)$","",after_stat(stratum))),
              size=1.5, fontface="bold") +
    scale_x_discrete(limits=c("From","To")) +
    scale_fill_manual(values=fill_vals, guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.01,0.03))) +
    theme_minimal(base_size=14) +
    theme(
      axis.title  = element_blank(),
      axis.text.x = element_text(face="bold", size=13),
      axis.text.y = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(5,15,5,5)
    )
}

build_backbone_plot <- function(bb, cluster_vec, pretty_vec, europe,
                                nrow_facet, fig_w, fig_h, out_path){
  p <- ggplot() +
    geom_sf(data=europe, fill="honeydew", color="black", linewidth=0.2) +
    geom_sf(data=bb$edges_ctx, color="grey75", linewidth=0.25, alpha=0.40) +
    geom_sf(data=bb$edges_top, aes(linewidth=edge_w),
            color="brown", alpha=0.95, lineend="round") +
    geom_sf(data=bb$nodes, aes(size=richness_mig, fill=betweenness),
            shape=21, color="black", stroke=0.25, alpha=0.95) +
    scale_linewidth_continuous(
      name   = "Corridor strength (relative within region)",
      breaks = c(0.7,1.2,1.7,2.2),
      labels = c("weak","moderate","strong","very strong"),
      range  = c(0.6,2.2)
    ) +
    scale_fill_distiller(
      name      = "Betweenness",
      palette   = "RdYlBu",
      direction = 1,
      limits    = c(0.01,0.06),
      breaks    = c(0.01,0.04,0.06),
      labels    = number_format(accuracy=0.01)
    ) +
    scale_size_continuous(
      name  = "Migration richness (unique genera)",
      range = c(1.5,5)
    ) +
    scale_x_continuous(breaks=c(-20,0,20,40),
                       labels=c("20°W","0°","20°E","40°E")) +
    coord_sf(xlim=c(-25,41), ylim=c(34,81), expand=FALSE) +
    facet_wrap(~GeoCluster, nrow=nrow_facet,
               labeller=labeller(GeoCluster=as_labeller(pretty_vec)),
               drop=FALSE) +
    guides(
      fill      = guide_colorbar(order=1, title.position="top", direction="horizontal",
                                 barwidth=unit(4,"cm"), barheight=unit(0.5,"cm")),
      size      = guide_legend(order=2, title.position="top"),
      linewidth = guide_legend(order=3, title.position="top",
                               override.aes=list(color="brown"))
    ) +
    theme_minimal(base_size=18) +
    theme(
      axis.text        = element_text(size=14),
      legend.title     = element_text(size=16, face="bold"),
      legend.text      = element_text(size=18),
      legend.position  = "top",
      legend.box       = "horizontal",
      panel.grid.major = element_line(color="grey92", linewidth=0.15),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face="bold", size=14),
      panel.spacing    = unit(0.8,"lines"),
      panel.border     = element_rect(colour="black", fill=NA, linewidth=1),
      plot.margin      = margin(t=2, r=10, b=10, l=10)
    )
  
  ggsave(out_path, p, width=fig_w, height=fig_h, dpi=300, units="in", bg="white")
  message("✓ Saved: ", out_path)
  p
}

out_10 <- file.path(out_dir, "season_migratory_raw")
out_4  <- file.path(out_dir, "4region_subset")
ensure_dir(out_10)
ensure_dir(out_4)

cluster_pretty_10_sankey <- c(
  "HC_Cluster_1_Alpine"        = "Central\nAlpine",
  "HC_Cluster_1_Atlantic"      = "Atlantic",
  "HC_Cluster_1_Continental"   = "Western\nContinental",
  "HC_Cluster_2_Alpine"        = "Eastern\nAlpine",
  "HC_Cluster_2_Continental"   = "Eastern\nContinental",
  "HC_Cluster_2_Mediterranean" = "Southeast\nMediterranean",
  "HC_Cluster_2_Pannonian"     = "Pannonian",
  "HC_Cluster_3_Alpine"        = "Scandinavian\nHighlands",
  "HC_Cluster_3_Boreal"        = "Boreal\nBaltic",
  "HC_Cluster_4_Mediterranean" = "Iberian"
)

cluster_pretty_4_sankey <- c(
  "HC_Cluster_1_Atlantic"    = "Atlantic",
  "HC_Cluster_1_Continental" = "Western\nContinental",
  "HC_Cluster_2_Continental" = "Eastern\nContinental",
  "HC_Cluster_3_Boreal"      = "Boreal\nBaltic"
)

message("• Plotting 10-cluster backbone …")
build_backbone_plot(
  bb         = bb10,
  cluster_vec = clusters_10,
  pretty_vec  = cluster_pretty_10,
  europe      = europe,
  nrow_facet  = 2,
  fig_w       = 26,
  fig_h       = 12,
  out_path    = file.path(out_10, "Anatidae_10cluster_Backbone_No_NS.png")
)

message("• Plotting 10-cluster Sankey …")
p_sank10 <- build_sankey(clusters_10, cluster_pretty_10_sankey, cluster_pal_10, M_SUM_NS)
ggsave(file.path(out_10, "Anatidae_10cluster_Sankey_No_NS.png"),
       p_sank10, width=12, height=8, dpi=300, units="in", bg="white")
message("✓ Saved: ", file.path(out_10, "Anatidae_10cluster_Sankey_No_NS.png"))

message("• Plotting 4-cluster backbone …")
build_backbone_plot(
  bb          = bb4,
  cluster_vec = clusters_4,
  pretty_vec  = cluster_pretty_4,
  europe      = europe,
  nrow_facet  = 1,
  fig_w       = 18,
  fig_h       = 8,
  out_path    = file.path(out_4, "Anatidae_4cluster_Backbone_No_NS.png")
)

message("• Plotting 4-cluster Sankey …")
p_sank4 <- build_sankey(clusters_4, cluster_pretty_4_sankey, cluster_pal_4, M_SUM_NS)
ggsave(file.path(out_4, "Anatidae_4cluster_Sankey_No_NS.png"),
       p_sank4, width=9, height=6, dpi=300, units="in", bg="white")
message("✓ Saved: ", file.path(out_4, "Anatidae_4cluster_Sankey_No_NS.png"))

message("• Plotting 4-cluster backbone (poster 2×2) …")
p_poster <- ggplot() +
  geom_sf(data=europe, fill="honeydew", color="black", linewidth=0.2) +
  geom_sf(data=bb4$edges_ctx, color="grey75", linewidth=0.25, alpha=0.40) +
  geom_sf(data=bb4$edges_top, aes(linewidth=edge_w),
          color="brown", alpha=0.95, lineend="round") +
  geom_sf(data=bb4$nodes, aes(size=richness_mig, fill=betweenness),
          shape=21, color="black", stroke=0.25, alpha=0.95) +
  scale_linewidth_continuous(
    name   = "Corridor strength",
    breaks = c(0.7,1.2,1.7,2.2),
    labels = c("weak","moderate","strong","very strong"),
    range  = c(0.6,2.2)
  ) +
  scale_fill_distiller(
    name      = "Betweenness",
    palette   = "RdYlBu",
    direction = 1,
    limits    = c(0.01,0.06),
    breaks    = c(0.01,0.04,0.06),
    labels    = number_format(accuracy=0.01)
  ) +
  scale_size_continuous(
    name  = "Migration richness (unique genera)",
    range = c(2,6)
  ) +
  scale_x_continuous(breaks=c(-20,0,20,40),
                     labels=c("20°W","0°","20°E","40°E")) +
  coord_sf(xlim=c(-25,41), ylim=c(34,81), expand=FALSE) +
  facet_wrap(~GeoCluster, nrow=2, ncol=2,
             labeller=labeller(GeoCluster=as_labeller(cluster_pretty_4)),
             drop=FALSE) +
  guides(
    fill      = guide_colorbar(order=1, title.position="top", direction="horizontal",
                               barwidth=unit(5,"cm"), barheight=unit(0.6,"cm")),
    size      = guide_legend(order=2, title.position="top"),
    linewidth = guide_legend(order=3, title.position="top",
                             override.aes=list(color="brown"))
  ) +
  theme_minimal(base_size=20) +
  theme(
    axis.text        = element_text(size=14),
    legend.title     = element_text(size=18, face="bold"),
    legend.text      = element_text(size=16),
    legend.position  = "top",
    legend.box       = "horizontal",
    panel.grid.major = element_line(color="grey92", linewidth=0.15),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face="bold", size=18),
    panel.spacing    = unit(1.2,"lines"),
    panel.border     = element_rect(colour="black", fill=NA, linewidth=1),
    plot.margin      = margin(t=5, r=15, b=10, l=10)
  )

ggsave(file.path(out_4, "Anatidae_4cluster_Backbone_poster_2x2_No_NS.png"),
       p_poster, width=12, height=12, dpi=300, units="in", bg="white")

