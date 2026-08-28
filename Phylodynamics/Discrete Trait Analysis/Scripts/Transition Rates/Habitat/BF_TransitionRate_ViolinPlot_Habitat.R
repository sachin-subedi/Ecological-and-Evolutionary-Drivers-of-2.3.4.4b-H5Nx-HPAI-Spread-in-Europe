options(scipen = 999)
library(dplyr)
library(stringr)
library(purrr)
library(coda)
library(readr)
library(ggplot2)
library(scales)

setwd(dir_logs)


burnin_pc       <- 50
reps            <- c("rep1", "rep2", "rep3")
sampling_labels <- c("equal", "proportional", "stratified")
BF_CUTOFF       <- 100    
LABEL_MODE      <- "columns"
POS_FROM  <- -0.58
POS_ARROW <- -0.52
POS_TO    <- -0.47
LEFT_MARGIN <- 60


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

file_manifest <- expand.grid(rep = reps, sampling = sampling_labels,
                             stringsAsFactors = FALSE) %>%
  mutate(
    log_file = paste0(rep, "_Habitat_rates_", sampling, "_combined.log")
  )

message("Files expected:")
walk(file_manifest$log_file, ~ message("  ", .x,
                                       if (file.exists(.x)) "  [found]" else "  [MISSING]"))

if (!any(file.exists(file_manifest$log_file)))
  stop("No log files found in:\n  ", getwd(),
       "\nSet dir_logs to the folder containing *_Habitat_rates_*_combined.log")

habitat_cols_solid <- c(
  "Coastal"        = "#0072B2",
  "Farm"           = "#CD5C5C",
  "Forest"         = "#BCBD22",
  "Grassland"      = "#26A69A",
  "Human_Modified" = "#56B4E9",
  "Marine"         = "#CC79A7",
  "Shrubland"      = "lightslategray",
  "Urban"          = "#E6842A",
  "Woodland"       = "#6B4C1B",
  "Rock"           = "#7570B3",
  "Wetland"        = "#006D2C"
)

lighten_palette <- function(cols, factor = 0.5) {
  out <- sapply(cols, function(color) {
    rgb_col <- col2rgb(color)
    new_col <- rgb_col + (255 - rgb_col) * factor
    rgb(new_col[1], new_col[2], new_col[3], maxColorValue = 255)
  })
  setNames(out, names(cols))
}
habitat_cols_light <- lighten_palette(habitat_cols_solid, factor = 0.5)

supported_pairs <- tibble::tribble(
  ~From,       ~To,
  "Wetland",   "Forest",
  "Wetland",   "Coastal",
  "Wetland",   "Farm",
  "Wetland",   "Grassland",
  "Grassland", "Coastal",
  "Grassland", "Farm",
  "Farm",      "Wetland",
  "Coastal",   "Forest",
  "Coastal",   "Grassland",
  "Coastal",   "Marine"
)

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


get_real_draws <- function(log_file, keep_pairs, burnin_pc = 50) {
  
  log_df <- read.delim(log_file, check.names = FALSE, sep = "\t")
  
  n_burn <- floor(burnin_pc / 100 * nrow(log_df))
  log_df <- slice(log_df, (n_burn + 1):n())
  
  rate_cols <- grep("^Habitat\\.rates\\.", names(log_df), value = TRUE)
  
  out <- map_dfr(rate_cols, function(rate_col) {
    
    from <- extract_second_last(rate_col)
    to   <- extract_last(rate_col)
    if (!any(keep_pairs$From == from & keep_pairs$To == to)) return(NULL)
    
    ind_col  <- gsub("rates", "indicators", rate_col, fixed = TRUE)
    rate_vec <- log_df[[rate_col]]
    ind_vec  <- log_df[[ind_col]]
    
    tibble(
      From      = from,
      To        = to,
      indicator = ind_vec,
      real_rate = rate_vec * ind_vec
    )
  })
  
  if (nrow(out) == 0)
    out <- tibble(From = character(), To = character(),
                  indicator = numeric(), real_rate = numeric())
  out
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
    out_file <- file.path(dir_out,
                          paste0("Habitat_bf_", rep_label, "_", samp_label, ".csv"))
    write_csv(result, out_file)
    message("  \u2705 Saved: ", basename(out_file))
    all_results[[length(all_results) + 1]] <- result
  }
}

if (length(all_results) == 0) stop("No results produced — check file names and paths.")

master_df <- bind_rows(all_results)
write_csv(master_df, file.path(dir_out, "Habitat_bf_ALL_reps_sampling.csv"))
message("\n\u2705 Master table saved  (", nrow(master_df), " rows, ",
        n_distinct(master_df$rep), " reps x ",
        n_distinct(master_df$sampling), " subsamples)")


n_legend <- list()

for (samp in sampling_labels) {
  
  message("\nPlotting: ", samp)
  
  ## (a) summary values
  df_bar <- master_df %>%
    filter(sampling == samp) %>%
    rename(From = from, To = to,
           Rate = median_rate, Lower = ci_lower, Upper = ci_upper) %>%
    semi_join(supported_pairs, by = c("From", "To")) %>%
    filter(bayes_factor >= BF_CUTOFF) %>%
    group_by(From, To) %>%
    summarise(Rate  = median(Rate),
              Lower = median(Lower),
              Upper = median(Upper), .groups = "drop") %>%
    arrange(desc(Rate))
  
  if (nrow(df_bar) == 0) {
    warning("No supported transitions for ", samp, " — skipped.")
    next
  }

  if (LABEL_MODE == "mono") {
    w_from <- max(nchar(df_bar$From))
    w_to   <- max(nchar(df_bar$To))
    make_lab <- function(f, t)
      paste0(str_pad(f, w_from, side = "left"), " \u2192 ",
             str_pad(t, w_to,   side = "right"))
  } else {
    ## "columns": row key is plain, text is drawn separately in aligned columns
    make_lab <- function(f, t) paste(f, "\u2192", t)
  }
  
  df_bar <- df_bar %>% mutate(Transition = make_lab(From, To))
  lev <- rev(df_bar$Transition)
  df_bar$Transition <- factor(df_bar$Transition, levels = lev)
  
  files <- paste0(reps, "_Habitat_rates_", samp, "_combined.log")
  ok    <- file.exists(files)
  if (!any(ok)) { warning("No logs for ", samp, " — skipped."); next }
  if (!all(ok)) warning("Missing: ", paste(files[!ok], collapse = ", "))
  
  df_v <- map_dfr(files[ok], ~ get_real_draws(.x, supported_pairs, burnin_pc)) %>%
    semi_join(df_bar %>% select(From, To), by = c("From", "To")) %>%
    mutate(Transition = factor(make_lab(From, To), levels = lev))
  
  if (nrow(df_v) == 0)
    stop("No draws extracted for ", samp,
         " — check that supported_pairs matches the log column names.")
  
  print(df_v %>%
          group_by(Transition) %>%
          summarise(pct_zero = round(100 * mean(indicator == 0), 1),
                    n = n(), .groups = "drop"))

  rng <- max(df_v$real_rate)
  
  p_violin <- ggplot() +
    geom_violin(data = df_v,
                aes(x = Transition, y = real_rate, fill = From),
                scale = "width", trim = TRUE, adjust = 1.0,
                colour = "grey30", linewidth = 0.4, alpha = 0.95) +
    geom_linerange(data = df_bar,
                   aes(x = Transition, ymin = Lower, ymax = Upper),
                   colour = "grey30", linewidth = 0.7) +
    geom_point(data = df_bar,
               aes(x = Transition, y = Rate),
               shape = 21, size = 2.6, stroke = 0.6,
               fill = "white", colour = "grey20")
  
  if (LABEL_MODE == "columns") {
    p_violin <- p_violin +
      geom_text(data = df_bar,
                aes(x = Transition, y = POS_FROM * rng, label = From),
                hjust = 1, size = 5.5, colour = "black") +
      geom_text(data = df_bar,
                aes(x = Transition, y = POS_ARROW * rng),
                label = "\u27F6", hjust = 0.5, size = 6.5, fontface = "bold", colour = "black") +
      geom_text(data = df_bar,
                aes(x = Transition, y = POS_TO * rng, label = To),
                hjust = 0, size = 5.5, colour = "black")
  }
  
  p_violin <- p_violin +
    coord_flip(clip = "off") +
    annotate("segment",
             x = 0.4, xend = length(lev) + 0.6,
             y = -0.6, yend = -0.6,
             colour = "black", linewidth = 0.8) +
    scale_fill_manual(values = habitat_cols_light,
                      limits = names(habitat_cols_light),
                      drop   = FALSE,
                      name   = "Source habitat") +
    scale_y_continuous(breaks = seq(0, ceiling(rng), by = 4),
                       labels = label_number(accuracy = 1),
                       expand = expansion(mult = c(0.38, 0.08))) +
    labs(x = "Habitat Transitions", y = "Transition rate") +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 16),
      axis.title.x = element_text(size = 18, face = "bold", hjust = 0.86,
                                  margin = margin(t = 12)),
      axis.title.y = element_text(size = 18, face = "bold", vjust = -10),
      axis.line.y  = element_blank()
      #axis.line.x  = element_blank()
    )
  
  if (LABEL_MODE == "columns") {
    p_violin <- p_violin +
      theme(axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin  = margin(t = 5, r = 12, b = 5, l = LEFT_MARGIN))
  } else {
    p_violin <- p_violin +
      theme(axis.text.y = element_text(size = 16, family = "mono", hjust = 1))
  }
  
  print(p_violin)
  
  out_name <- file.path(dir_out,
                        paste0("Habitat_Rates_violin_", samp, "_light0.5_CI.png"))
  ggsave(out_name, p_violin, width = 8, height = 10, bg = "white", dpi = 300)
  ggsave(sub("\\.png$", ".pdf", out_name), p_violin,
         width = 8, height = 10, device = cairo_pdf)
  message("  \u2705 saved: ", basename(out_name))
  
  n_legend[[samp]] <- df_v %>%
    count(Transition, name = "n_posterior_samples") %>%
    mutate(sampling = samp)
}


n_table <- bind_rows(n_legend)
write_csv(n_table, file.path(dir_out, "Fig3A_violin_n_per_transition.csv"))
print(n_table)


master_df %>%
  rename(From = from, To = to) %>%
  semi_join(supported_pairs, by = c("From", "To")) %>%
  filter(bayes_factor >= BF_CUTOFF) %>%
  write_csv(file.path(dir_out, "SourceData_Fig3A.csv"))

Try