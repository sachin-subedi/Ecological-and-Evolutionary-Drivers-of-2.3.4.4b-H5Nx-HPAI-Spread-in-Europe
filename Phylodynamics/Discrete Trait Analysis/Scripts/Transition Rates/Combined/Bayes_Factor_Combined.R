setwd("results/rates/Combined/")

library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(readr)

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
  
  log_df <- read.delim(log_file, check.names = FALSE, sep = "\t")
  
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- dplyr::slice(log_df, (n_burn + 1):dplyr::n())
  
  rate_cols <- grep("^HG\\.rates\\.", names(log_df), value = TRUE)
  
  k <- find_k(length(rate_cols))
  if (is.null(k)) stop("Could not infer k from number of rate columns in ", log_file)
  
  q_prior <- (log(2) + k - 1) / (k * (k - 2) / 2)
  
  summaries <- map_dfr(rate_cols, function(rate_col) {
    
    ind_col  <- gsub("rates", "indicators", rate_col, fixed = TRUE)
    rate_vec <- log_df[[rate_col]]
    ind_vec  <- log_df[[ind_col]]
    real_vec <- rate_vec * ind_vec
    
    tibble(
      from           = extract_second_last(rate_col),
      to             = extract_last(rate_col),
      mean_indicator = mean(ind_vec),
      mean_rate      = mean(rate_vec),
      mean_real_rate = mean(real_vec),
      median_rate    = median(real_vec),
      hpd_lower      = HPDinterval(as.mcmc(real_vec), prob = 0.95)[1, "lower"],
      hpd_upper      = HPDinterval(as.mcmc(real_vec), prob = 0.95)[1, "upper"],
      bayes_factor   = (mean(ind_vec) * (1 - q_prior)) /
        ((1 - mean(ind_vec)) * q_prior)
    )
  })
  
  summaries
}

datasets <- list(
  list(log_file = "HG_rates_equal_combined.log",        out_file = "HG_bf_equal.csv"),
  list(log_file = "HG_rates_proportional_combined.log", out_file = "HG_bf_proportional.csv"),
  list(log_file = "HG_rates_stratified_combined.log",   out_file = "HG_bf_stratified.csv")
)

for (ds in datasets) {
  message("Processing: ", ds$log_file)
  summaries <- process_log_file(ds$log_file, burnin_pc = burnin_pc)
  write_csv(summaries, ds$out_file)
  message("✓ ", ds$out_file, " saved.")
}
