library(dplyr)
library(readr)

setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/reproduce/rates/Combined")
df <- read_csv("HG_bf_Subsample1.csv")

df_rounded <- df %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(is.infinite(.x), .x, round(.x, 2))
    )
  )

write_csv(df_rounded, "HG_bf_Subsample1_dec.csv")
