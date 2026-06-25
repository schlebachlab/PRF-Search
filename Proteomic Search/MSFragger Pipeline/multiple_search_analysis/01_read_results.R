library(tidyverse)

#for loop to pull in classified peptides from the different searches
search_conditions <- list("main_search", "downstream_search", "upstream_search")

for (condition in search_conditions) {
  directory <- paste0("../second_pass_", condition)
  file_path <- paste0(directory, "/analysis/output_tables/classified_peptides_data.tsv")
  loop_tbl <- read_tsv(file_path, col_types = cols("prev_aa" = col_character(), "charges" = col_character()))
  assign(paste0(condition, "_classified_peptides"), loop_tbl, envir = .GlobalEnv)
  rm(directory, file_path, loop_tbl)
}
