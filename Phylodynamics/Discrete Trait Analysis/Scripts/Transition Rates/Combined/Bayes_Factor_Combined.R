library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(readr)
library(data.table)
library(writexl)

find_k <- function(m) {
  for (n in 1:m) {
    if (m == n * (n - 1)) return(n)
  }
  NULL
}

extract_second_last <- function(x) {
  parts <- str_split(x, "\\.", simplify = TRUE)
  parts[, ncol(parts) - 1]
}

extract_last <- function(x) {
  parts <- str_split(x, "\\.", simplify = TRUE)
  parts[, ncol(parts)]
}

burnin_pc <- 50

process_log_file <- function(log_file, burnin_pc = 50) {
  
  header <- names(fread(log_file, nrows = 0, header = TRUE, sep = "\t"))
  
  keep_cols <- grep("^HG\\.(rates|indicators)\\.", header, value = TRUE)
  if (length(keep_cols) == 0)
    stop("No HG.rates/HG.indicators columns found in ", log_file)
  
  log_df <- fread(log_file, sep = "\t", header = TRUE,
                  select = keep_cols, data.table = FALSE)
  
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- log_df[(n_burn + 1):nrow(log_df), , drop = FALSE]
  
  rate_cols <- grep("^HG\\.rates\\.", names(log_df), value = TRUE)
  
  k <- find_k(length(rate_cols))
  if (is.null(k)) stop("Could not infer k from number of rate columns in ", log_file)
  
  q_prior <- (log(2) + k - 1) / (k * (k - 2) / 2)
  
  summaries <- map_dfr(rate_cols, function(rate_col) {
    
    ind_col  <- gsub("rates", "indicators", rate_col, fixed = TRUE)
    rate_vec <- log_df[[rate_col]]
    ind_vec  <- log_df[[ind_col]]
    real_vec <- rate_vec * ind_vec
    
    hpd <- HPDinterval(as.mcmc(real_vec), prob = 0.95)
    
    tibble(
      from           = extract_second_last(rate_col),
      to             = extract_last(rate_col),
      mean_indicator = mean(ind_vec),
      mean_rate      = mean(rate_vec),
      mean_real_rate = mean(real_vec),
      median_rate    = median(real_vec),
      hpd_lower      = HPDinterval(as.mcmc(real_vec), prob = 0.95)[1, "lower"],
      hpd_upper      = HPDinterval(as.mcmc(real_vec), prob = 0.95)[1, "upper"],
      ci_lower       = unname(quantile(real_vec, 0.025)),   # equal-tailed 95% credible interval
      ci_upper       = unname(quantile(real_vec, 0.975)),   # equal-tailed 95% credible interval
      bayes_factor   = (mean(ind_vec) * (1 - q_prior)) /
        ((1 - mean(ind_vec)) * q_prior)
    )
  })
  
  summaries <- summaries %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
  
  summaries
}

reps   <- c("rep1", "rep2", "rep3")
models <- c("equal", "proportional", "stratified")

datasets <- list()
for (rep in reps) {
  for (model in models) {
    datasets[[length(datasets) + 1]] <- list(
      rep      = rep,
      model    = model,
      log_file = sprintf("%s_HG_rates_%s_combined.log", rep, model),
      out_file = sprintf("%s_HG_bf_%s.csv", rep, model)
    )
  }
}

all_summaries <- list()

for (ds in datasets) {
  message("Processing: ", ds$log_file)
  summaries <- process_log_file(ds$log_file, burnin_pc = burnin_pc)
  
  write_csv(summaries, ds$out_file)
  message("✓ ", ds$out_file, " saved.")
  
  all_summaries[[length(all_summaries) + 1]] <- summaries %>%
    mutate(subsampling = ds$model, rep = ds$rep, .before = 1)
}

combined_df <- bind_rows(all_summaries)

write_xlsx(combined_df, "HG_bf_ALL_combined.xlsx")

region_labels <- c(
  "HC1_Alp" = "Central Alpine",
  "HC1_Atl" = "Atlantic",
  "HC1_Con" = "Western Continental",
  "HC2_Alp" = "Eastern Alpine",
  "HC2_Con" = "Eastern Continental",
  "HC2_Med" = "Southeast Mediterranean",
  "HC2_Pan" = "Pannonian",
  "HC3_Alp" = "Scandinavian Highlands",
  "HC3_Bor" = "Boreal Baltic",
  "HC4_Med" = "Iberian"
)

hab_map <- c(
  CM = "Coastal",
  FA = "Farm",
  FO = "Forest",
  GW = "Grassland",
  WT = "Wetland",
  UB = "Urban"
)

recode_label <- function(x) {
  x_chr <- as.character(x)
  sapply(x_chr, function(val) {
    parts    <- strsplit(val, "_")[[1]]
    hab_code <- tail(parts, 1)                              
    geo_code <- paste(parts[-length(parts)], collapse = "_") 
    geo_full <- region_labels[geo_code]
    hab_full <- hab_map[hab_code]
    if (!is.na(geo_full) && !is.na(hab_full)) {
      paste(geo_full, hab_full)               
    } else {
      val
    }
  }, USE.NAMES = FALSE)
}

combined_df <- combined_df %>%
  mutate(
    from = recode_label(from),
    to   = recode_label(to)
  )

write_xlsx(combined_df, "HG_bf_ALL_combined_decoded.xlsx")

