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

burnin_pc       <- 50
reps            <- c("rep1", "rep2", "rep3")
sampling_labels <- c("equal", "proportional", "stratified")

file_manifest <- expand.grid(rep = reps, sampling = sampling_labels,
                             stringsAsFactors = FALSE) %>%
  mutate(
    log_file = paste0(rep, "_Habitat_rates_", sampling, "_combined.log")
  )

message("Files expected:")
walk(file_manifest$log_file, ~ message("  ", .x,
                                       if (file.exists(.x)) "  [found]" else "  [MISSING]"))

process_log_file <- function(log_file, burnin_pc = 50) {
  
  log_df <- read.delim(log_file, check.names = FALSE, sep = "\t")
  
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- slice(log_df, (n_burn + 1):n())
  
  rate_cols <- grep("^Habitat\\.rates\\.", names(log_df), value = TRUE)
  
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
      ci_lower       = unname(quantile(real_vec, 0.025)),
      ci_upper       = unname(quantile(real_vec, 0.975)),
      bayes_factor   = (mean(ind_vec) * (1 - q_prior)) /
        ((1 - mean(ind_vec)) * q_prior)
    )
  })
  
  summaries %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

all_results <- list()

for (i in seq_len(nrow(file_manifest))) {
  
  log_file   <- file_manifest$log_file[i]
  rep_label  <- file_manifest$rep[i]
  samp_label <- file_manifest$sampling[i]
  
  if (!file.exists(log_file)) {
    warning("File not found, skipping: ", log_file)
    next
  }
  
  message("Processing: ", log_file)
  
  result <- tryCatch(
    process_log_file(log_file, burnin_pc) %>%
      mutate(rep = rep_label, sampling = samp_label, .before = 1),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  
  if (!is.null(result)) {
    # per-file CSV
    out_file <- paste0("Habitat_bf_", rep_label, "_", samp_label, ".csv")
    write_csv(result, out_file)
    message("  \u2705 Saved: ", out_file)
    
    all_results[[length(all_results) + 1]] <- result
  }
}

if (length(all_results) > 0) {
  master_df <- bind_rows(all_results)
  write_csv(master_df, "Habitat_bf_ALL_reps_sampling.csv")
  message("\n\u2705 Master table saved: Habitat_bf_ALL_reps_sampling.csv",
          "  (", nrow(master_df), " rows, ",
          n_distinct(master_df$rep), " reps x ",
          n_distinct(master_df$sampling), " subsamples)")
} else {
  message("No results produced — check file names and paths.")
}
