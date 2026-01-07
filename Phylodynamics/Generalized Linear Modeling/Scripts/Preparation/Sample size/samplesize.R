library(dplyr)
library(readr)
library(lubridate)

setwd("samplesize/")
subsample_data <- read.delim("subsampled_data.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(
    Collection_Date = as.Date(Collection_Date),
    Year = year(Collection_Date)
  )

GeoCluster_sample_size <- subsample_data %>%
  group_by(GeoCluster) %>%
  summarise(sample_size = n(), .groups = "drop") %>%
  arrange(desc(sample_size))

write_tsv(GeoCluster_sample_size, "GeoCluster_sample_size.tsv")

cat("✅ File 'sample_size.tsv' has been saved.\n")
