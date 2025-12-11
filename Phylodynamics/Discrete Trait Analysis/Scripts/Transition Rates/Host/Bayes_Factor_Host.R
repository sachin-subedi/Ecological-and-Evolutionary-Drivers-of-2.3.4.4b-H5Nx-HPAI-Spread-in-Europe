##############################################################################
# 0. Libraries ---------------------------------------------------------------
##############################################################################
library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(readr)

##############################################################################
# 1. Helper utilities --------------------------------------------------------
##############################################################################
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

##############################################################################
# 2. File paths & settings ---------------------------------------------------
##############################################################################
workdir   <- "~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/rates/Host"
setwd(workdir)

burnin_pc <- 50   # % burn-in to discard

log_files <- c(
  "Host_rates_Subsample1_combined.log",
  "Host_rates_Subsample2_combined.log",
  "Host_rates_Subsample3_combined.log"
)

##############################################################################
# 3. Function to process ONE log file ---------------------------------------
##############################################################################
process_log_file <- function(log_file, burnin_pc = 50) {
  
  # read log
  log_df <- read.delim(log_file, check.names = FALSE, sep = "\t")
  
  # apply burn-in
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- dplyr::slice(log_df, (n_burn + 1):dplyr::n())
  
  # identify rate columns
  rate_cols <- grep("^Host\\.rates\\.", names(log_df), value = TRUE)
  
  # infer k and prior inclusion probability
  k <- find_k(length(rate_cols))
  if (is.null(k)) stop("Could not infer k from number of rate columns in ", log_file)
  
  q_prior <- (log(2) + k - 1) / (k * (k - 2) / 2)
  
  # per-transition summaries
  summaries <- map_dfr(rate_cols, function(rate_col){
    
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

##############################################################################
# 4. Loop over subsample logs & save outputs --------------------------------
##############################################################################
for (log_file in log_files) {
  
  # e.g. "Host_rates_Subsample1_combined.log" -> "Subsample1"
  subs_label <- gsub("Host_rates_(Subsample[0-9]+)_combined\\.log", "\\1", log_file)
  
  message("Processing: ", log_file, "  (", subs_label, ")")
  
  summaries <- process_log_file(log_file, burnin_pc = burnin_pc)
  
  out_file <- paste0("Host_bf_", subs_label, ".csv")
  write_csv(summaries, out_file)
  
  message("✓ ", out_file, " saved.")
}

##Combined
all_summaries <- map2_dfr(
  log_files,
  c("Subsample1", "Subsample2", "Subsample3"),
  ~ process_log_file(.x, burnin_pc) %>% mutate(subsample = .y)
)

write_csv(all_summaries, "Host_bf_all_subsamples.csv")
