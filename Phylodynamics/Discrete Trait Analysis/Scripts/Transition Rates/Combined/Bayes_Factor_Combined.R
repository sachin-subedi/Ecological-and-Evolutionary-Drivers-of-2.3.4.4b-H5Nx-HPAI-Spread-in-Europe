library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(readr)

find_k <- function(m){
  for (n in 1:m) {
    if (m == n * (n - 1)) return(n)
  }
  NULL
}

extract_second_last <- function(x){
  parts <- str_split(x, "\\.", simplify = TRUE)
  parts[, ncol(parts) - 1]
}

extract_last <- function(x){
  parts <- str_split(x, "\\.", simplify = TRUE)
  parts[, ncol(parts)]
}

rename_hab_code <- function(x) {
  y <- as.character(x)
  y <- gsub("^CMW(\\d)$",    "Coastal\\1",   y)
  y <- gsub("^FD(\\d)$",     "Farm\\1",      y)
  y <- gsub("^GW(\\d)$",     "Grassland\\1", y)
  y <- gsub("^RSFWW(\\d)$",  "Forest\\1",    y)
  y <- gsub("^WW(\\d)$",     "Wetland\\1",   y)
  y
}

workdir   <- "~/Combined"
setwd(workdir)

burnin_pc <- 50

log_files <- c(
  "HG_rates_Subsample1_combined.log",
  "HG_rates_Subsample2_combined.log",
  "HG_rates_Subsample3_combined.log"
)

process_log_file <- function(log_file, burnin_pc = 50) {
  
  log_df <- read.delim(log_file, check.names = FALSE, sep = "\t")
  
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- dplyr::slice(log_df, (n_burn + 1):dplyr::n())
  
  rate_cols <- grep("^HG\\.rates\\.", names(log_df), value = TRUE)
  
  k <- find_k(length(rate_cols))
  if (is.null(k)) stop("Could not infer k from number of rate columns in ", log_file)
  
  q_prior <- (log(2) + k - 1) / (k * (k - 2) / 2)
  
  summaries <- map_dfr(rate_cols, function(rate_col){
    
    ind_col  <- gsub("rates", "indicators", rate_col, fixed = TRUE)
    rate_vec <- log_df[[rate_col]]
    ind_vec  <- log_df[[ind_col]]
    real_vec <- rate_vec * ind_vec
    
    from_code <- extract_second_last(rate_col)
    to_code   <- extract_last(rate_col)
    
    tibble(
      from           = rename_hab_code(from_code),
      to             = rename_hab_code(to_code),
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

for (log_file in log_files) {
  
  subs_label <- gsub("HG_rates_(Subsample[0-9]+)_combined\\.log", "\\1", log_file)
  
  message("Processing: ", log_file, "  (", subs_label, ")")
  
  summaries <- process_log_file(log_file, burnin_pc = burnin_pc)
  
  out_file <- paste0("HG_bf_", subs_label, ".csv")
  write_csv(summaries, out_file)
  
  message("✓ ", out_file, " saved.")
}

all_summaries <- map2_dfr(
  log_files,
  c("Subsample1", "Subsample2", "Subsample3"),
  ~ process_log_file(.x, burnin_pc) %>% mutate(subsample = .y)
)

write_csv(all_summaries, "HG_bf_all_subsamples.csv")

