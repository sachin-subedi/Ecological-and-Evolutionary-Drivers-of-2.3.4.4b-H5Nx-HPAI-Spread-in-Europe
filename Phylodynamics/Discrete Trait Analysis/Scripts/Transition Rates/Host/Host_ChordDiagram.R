library(circlize)
library(tidytree)
library(RColorBrewer)

options(scipen = 999)

setwd("~/Host")
custom_colors <- c(
  DomesticMammal = "#6A3D9A",
  DomesticBird   = "#CD5C5C",
  WildMammal     = "#BCBD22",
  Human          = "#1B9E77",
  WildBird       = "#3C5488"
)

bayes_colors <- c(
  ">100"   = "grey0",
  "31-100" = "grey60",
  "11-30"  = "grey80",
  "3-10"   = "grey90",
  "<3"     = "grey95"
)

bayes_labels_clean <- c(
  "Decisive (>100)",
  "Very Strong (31–100)",
  "Strong (11–30)",
  "Substantial (3–10)"
)
host_files <- c(
  "Host_bf_Subsample1.csv",
  "Host_bf_Subsample2.csv",
  "Host_bf_Subsample3.csv"
)
sub_labels <- c("Subsample1", "Subsample2", "Subsample3")

for (i in seq_along(host_files)) {
  in_file <- host_files[i]
  tag     <- sub_labels[i]
  message("Processing ", in_file, " (", tag, ")")
  
  df <- read.csv(in_file, header = TRUE)
  
  df$FROM <- gsub("_", "", as.character(df$from))
  df$TO   <- gsub("_", "", as.character(df$to))
  df$FROM <- gsub("\\s+", "", df$FROM)
  df$TO   <- gsub("\\s+", "", df$TO)
  
  df$FROM <- paste0("src-", df$FROM)
  df$TO   <- paste0("snk-", df$TO)
  
  if ("bf" %in% names(df)) {
    BF_raw <- df$bf
  } else if ("bayes_factor" %in% names(df)) {
    BF_raw <- df$bayes_factor
  } else {
    stop("No 'bf' or 'bayes_factor' column found in ", in_file)
  }
  
  df$BF <- as.numeric(as.character(BF_raw))
  
  rate_col <- NULL
  if ("posterior.mean.of.Host.rates" %in% names(df)) {
    rate_col <- "posterior.mean.of.Host.rates"
  } else if ("mean_rate" %in% names(df)) {
    rate_col <- "mean_rate"
  } else {
    stop("No rate column ('posterior.mean.of.Host.rates' or 'mean_rate') in ", in_file)
  }
  
  df$BFcategory <- with(df,
                        ifelse(BF >= 100, ">100",
                               ifelse(BF >= 31, "31-100",
                                      ifelse(BF >= 11, "11-30",
                                             ifelse(BF >= 3, "3-10", "<3")))))
  
  df$BFcategory <- factor(df$BFcategory,
                          levels = c(">100", "31-100", "11-30", "3-10", "<3"))
  
  df$BAYES_FACTOR <- bayes_colors[df$BFcategory]
  
  Host_names  <- gsub("src-|snk-", "", unique(c(df$FROM, df$TO)))
  Host_colors <- custom_colors[Host_names]
  names(Host_colors) <- Host_names
  
  grid.col <- setNames(
    Host_colors[gsub("src-|snk-", "", unique(c(df$FROM, df$TO)))],
    unique(c(df$FROM, df$TO))
  )
  
  unique_sectors <- unique(c(df$FROM, df$TO))
  gap.after <- rep(1, length(unique_sectors))
  gap.after[length(unique(df$FROM))] <- 30
  gap.after[length(unique_sectors)]  <- 30
  
  circos.clear()
  circos.par(start.degree = 90, gap.after = gap.after)
  
  reversed_Hosts <- rev(unique_sectors)
  
  png(
    file = paste0(
      "Host_ChordDiagram_",
      tag, ".png"
    ),
    width = 1800,
    height = 1200,
    res = 300,
    bg = "transparent"
  )
  
  par(
    mfrow  = c(1, 1),
    cex    = 0.8,
    family = "serif",
    mar    = c(1.5, 4, 2.5, 4)
  )
  
  chordDiagram(
    df[, c("FROM", "TO", rate_col)],
    order             = reversed_Hosts,
    grid.col          = grid.col,
    col               = df$BAYES_FACTOR,
    transparency      = 0,
    directional       = 1,
    direction.type    = "arrows",
    link.arr.type     = "big.arrow",
    annotationTrack   = "grid",
    preAllocateTracks = 1
  )
  
  lines(x = c(0, 0), y = c(-0.8, 1), lty = 2, col = "black", lwd = 1.2)
  
  text(-0.03, 0.95, "Source", cex = 1.3, font = 2, adj = 1)
  text( 0.03, 0.95, "Sink",   cex = 1.3, font = 2, adj = 0)
  
  legend_cols <- bayes_colors[c(">100", "31-100", "11-30", "3-10")]
  
  legend(
    "top",
    inset      = c(0, -0.07),
    legend     = bayes_labels_clean,
    col        = legend_cols,
    pch        = 15,
    cex        = 0.85,
    bty        = "n",
    xpd        = TRUE,
    ncol       = 4
  )

  dev.off()
}

