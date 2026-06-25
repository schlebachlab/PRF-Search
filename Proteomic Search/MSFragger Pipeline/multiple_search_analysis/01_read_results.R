# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

#for loop to pull in classified peptides from the different searches
search_conditions <- list("main_search", "downstream_search", "upstream_search")

for (condition in search_conditions) {
  directory <- paste0("../second_pass_", condition)
  file_path <- paste0(directory, "/analysis/output_tables/classified_peptides_data.tsv")
  loop_tbl <- read_tsv(file_path, col_types = cols("prev_aa" = col_character(),
                                                   "next_aa" = col_character(),
                                                   "frame_transition" = col_character(),
                                                   "charges" = col_character())) |> 
    rename(enz_act = source)
  assign(paste0(condition, "_classified_peptides"), loop_tbl, envir = .GlobalEnv)
  rm(directory, file_path, loop_tbl, condition)
}

all_search_classified_peptides <- bind_rows(upstream = upstream_search_classified_peptides, 
                                            prf_search = main_search_classified_peptides, 
                                            downstream = downstream_search_classified_peptides,
                                            .id = "source")

# pull data from the first pass search analysis
first_pass_chunks <- list() # create empty list to collect the chunks of first pass results

for (chunk in 1:3) {
  file_path <- paste0("../first_pass_search/analysis/first_pass_peptide_and_genes_chunk", chunk, ".tsv")
  first_pass_chunks[[chunk]] <- read_tsv(
    file_path, 
    col_types = cols("prev_aa" = col_character(), 
                     "next_aa" = col_character(), 
                     "charges" = col_character()
                     )
    )
  rm(file_path)
} # for loop to read each chunk

first_pass_peptides_and_genes <- bind_rows(first_pass_chunks)
