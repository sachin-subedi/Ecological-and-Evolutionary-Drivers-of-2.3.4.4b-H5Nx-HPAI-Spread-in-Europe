library(dplyr)
library(readr)
library(stringr)

var_name_mapping <- c(
  "1" = "Live Poultry (NetWeight)",
  "2" = "Sample Size Origin",
  "3" = "Sample Size Destination",
  "4" = "Wetland Areas Origin",
  "5" = "Wetland Areas Destination",
  "6" = "Temperature Origin",
  "7" = "Temperature Destination",
  "8" = "Anatidae Migration"
)

calc_min_interval <- function(x, alpha = 0.05) {
  x <- sort(x)
  n <- length(x)
  cred_mass <- 1 - alpha
  interval_idx_inc <- floor(cred_mass * n)
  n_intervals <- n - interval_idx_inc
  widths <- x[(interval_idx_inc + 1):n] - x[1:n_intervals]
  if (length(widths) == 0) stop("Too few elements for interval calculation")
  min_idx <- which.min(widths)
  c(x[min_idx], x[min_idx + interval_idx_inc])
}

BayesFactor <- function(indicators, n) {
  prior <- 1 - exp(log(0.5) / n)
  posterior <- mean(indicators, na.rm = TRUE)
  if (posterior == 1) posterior <- (length(indicators) - 1) / length(indicators)
  (posterior / (1 - posterior)) / (prior / (1 - prior))
}

GLMSummarize <- function(dirname) {
  setwd(dirname)
  files <- list.files(pattern = "\\.log$")
  results <- list()
  
  for (filename in files) {
    if (!grepl("glm", filename)) next
    message("Processing: ", filename)
    data <- read_tsv(filename, show_col_types = FALSE)
    model <- tools::file_path_sans_ext(filename)
    varNum <- as.integer(str_extract(names(data)[ncol(data)], "\\d+"))
    
    time_cols <- grep("Times", names(data), value = TRUE)
    for (time_col in time_cols) {
      varID <- str_extract(time_col, "\\d+")
      ind_cols <- grep(paste0("coefIndicators", varID), names(data), value = TRUE)
      if (length(ind_cols) == 0) next
      ind_col <- ind_cols[1]
      trait <- str_extract(ind_col, "^[^.]+")
      indicators <- data[[ind_col]]
      
      if (max(indicators, na.rm = TRUE) == 0) {
        lower_hpd <- 0; upper_hpd <- 0; median_val <- 0
      } else {
        vals <- data[[time_col]]
        vals[vals == 0] <- NA
        vals <- na.omit(vals)
        interval <- calc_min_interval(vals)
        lower_hpd <- interval[1]; upper_hpd <- interval[2]
        median_val <- median(vals)
      }
      
      pp <- mean(indicators, na.rm = TRUE)
      bf <- BayesFactor(indicators, varNum)
      
      desc <- if (varID %in% names(var_name_mapping)) var_name_mapping[[varID]] else varID
      
      results[[length(results) + 1]] <- data.frame(
        Model = model,
        Trait = trait,
        Variable = desc,
        CoefMedian = median_val,
        CoefLowerHPD = lower_hpd,
        CoefUpperHPD = upper_hpd,
        pp = pp,
        BF = bf,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results) > 0) {
    out_df <- bind_rows(results)
    write.table(out_df,
                file = file.path(dirname, "Subsample1_Set1a_GLMSummary.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    message("✅ Summary written to SoCluster_GLMSummary.txt")
  } else {
    warning("No results to save.")
  }
}

dirname <- "/Set1a/"
GLMSummarize(dirname)

