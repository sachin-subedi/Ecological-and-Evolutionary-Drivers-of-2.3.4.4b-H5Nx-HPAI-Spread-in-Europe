setwd("~/Library/CloudStorage/OneDrive-UniversityofGeorgia/EU_H5/may/data/Again/after_draft/pipeline/subsampled3/emptrees/")
subsample_trees <- function(input_file, output_file, n_trees = 500) {
  lines <- readLines(input_file)
  
  tree_indices <- which(grepl("tree STATE", lines))
  end_index <- which(grepl("End;", lines))[1]
  start_index <- tree_indices[1]
  
  before_tree_state <- lines[1:(start_index - 1)]
  tree_lines <- lines[tree_indices]
  
  n_sample <- min(n_trees, length(tree_lines))
  sampled_trees <- sample(tree_lines, n_sample)
  
  writeLines(c(before_tree_state, sampled_trees, "End;"), output_file)
}

count_trees <- function(tree_file) {
  lines <- readLines(tree_file)
  sum(grepl("tree STATE", lines))
}

input_path <- "empirical_subsampled3.trees"
output_path <- "empirical_subsampled3_500.trees"

subsample_trees(input_path, output_path, n_trees = 500)
# Verify tree count
num_trees <- count_trees(output_path)
cat("✅ Number of trees in output:", num_trees, "\n")
