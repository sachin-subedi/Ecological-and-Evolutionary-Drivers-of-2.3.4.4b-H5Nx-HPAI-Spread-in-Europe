setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/subsampled3/combined/")

data1 <- read.delim("subsampled_data.tsv", sep="\t", stringsAsFactors=FALSE)
colnames(data1)
library(dplyr)
library(stringr)
library(readr)

parse_cluster_num <- function(x) {
  tok <- x %>%
    str_replace("^\\s*GeoCluster[_\\s-]*", "") %>%
    str_trim()
  
  case_when(
    tok %in% c("1","2","3","4","5") ~ tok,
    str_to_lower(tok) == "one"   ~ "1",
    str_to_lower(tok) == "two"   ~ "2",
    str_to_lower(tok) == "three" ~ "3",
    str_to_lower(tok) == "four"  ~ "4",
    str_to_lower(tok) == "five"  ~ "5",
    TRUE ~ NA_character_
  )
}

data1 <- data1 %>%
  mutate(
    habitat_raw = Habitat %>% str_trim(),
    
    habitat_code = case_when(
      str_detect(habitat_raw, regex("^Wetland(s)?$", ignore_case = TRUE)) ~ "WW",
      str_detect(habitat_raw, regex("^Farm$",        ignore_case = TRUE)) ~ "FD",
      str_detect(habitat_raw, regex("^Grassland$",   ignore_case = TRUE)) ~ "GW",
      str_detect(habitat_raw, regex("^(Human_Modified|Rock|Shrubland|Forest|Woodland)$", ignore_case = TRUE)) ~ "RSFWW",
      str_detect(habitat_raw, regex("^Urban$",       ignore_case = TRUE)) ~ "UH",
      str_detect(habitat_raw, regex("^(Coastal|Marine)$", ignore_case = TRUE)) ~ "CMW",
      TRUE ~ "O"  # Others
    ),
    
    cluster_num = parse_cluster_num(GeoCluster),
    
    HG_compact = if_else(!is.na(cluster_num), paste0(habitat_code, cluster_num), NA_character_)
  ) %>%
  select(-habitat_raw)

write_tsv(data1, "combined_compact_subsampled_data3.tsv")

## Fasta file
library(readr)
library(stringr)
library(dplyr)

base_dir <- "~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/subsampled3/combined/"
tsv_path  <- file.path(base_dir, "combined_compact_subsampled_data3.tsv")
fasta_out <- file.path(base_dir, "combined_compact_subsampled_data3.fasta")

subsampled_data1 <- read.delim(
  tsv_path, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)

write_fasta <- function(df, out_path,
                        header_fields = c("Isolate_Id","Isolate_Name","Habitat","Type",
                                          "Country","GeoCluster","HG_compact","Collection_Date_Final"),
                        seq_field = "Sequence",
                        line_width = 60,
                        skip_empty_seq = TRUE) {
  
  missing_cols <- setdiff(header_fields, names(df))
  if (length(missing_cols) > 0) {
    df[missing_cols] <- NA_character_
  }
  if (!seq_field %in% names(df)) {
    stop(sprintf("Sequence field '%s' is missing.", seq_field))
  }
  
  con <- file(out_path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  n_written <- 0L
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    
    vals <- vapply(header_fields, function(nm) {
      v <- row[[nm]]
      if (is.null(v) || is.na(v) || (is.character(v) && nchar(trimws(v)) == 0)) "NA" else as.character(v)
    }, FUN.VALUE = character(1))
    
    header <- paste0(">", paste(vals, collapse = "|"))
    
    seq_raw <- row[[seq_field]]
    if (is.null(seq_raw) || is.na(seq_raw)) seq_raw <- ""
    seq_str <- gsub("\\s+", "", as.character(seq_raw))
    
    if (skip_empty_seq && nchar(seq_str) == 0) next
    
    if (nchar(seq_str) > 0 && is.finite(line_width) && line_width > 0) {
      seq_lines <- strwrap(seq_str, width = line_width)
    } else {
      seq_lines <- seq_str
    }
    
    writeLines(header, con)
    writeLines(seq_lines, con)
    n_written <- n_written + 1L
  }
  
  message(sprintf("Wrote %d FASTA records to: %s", n_written, out_path))
}

write_fasta(subsampled_data1, fasta_out)
cat("FASTA file saved as:", fasta_out, "\n")


