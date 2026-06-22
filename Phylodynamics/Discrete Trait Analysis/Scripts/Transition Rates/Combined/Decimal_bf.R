setwd("combined_traits/results/rates/sig")

library(dplyr)
library(readr)

cols <- c("from","to","mean_indicator","mean_rate","median_rate",
          "hpd_lower","hpd_upper","bayes_factor","subsample")

process_file <- function(infile, outbase) {
  df <- read_csv(infile, show_col_types = FALSE)
  if (!"subsample" %in% names(df) && "sampling" %in% names(df)) {
    df <- df %>% rename(subsample = sampling)
  }
  df <- df %>%
    select(all_of(cols)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write_csv(df, paste0(outbase, "_sig.csv"))
}

process_file("GeoCluster_bf_combined.csv", "GeoCluster_bf_combined")
process_file("Habitat_bf_all_sampling.csv", "Habitat_bf_all_sampling")
process_file("HG_bf_combined.csv", "HG_bf_combined")

