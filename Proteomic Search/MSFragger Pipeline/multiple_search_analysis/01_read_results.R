# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

#for loop to pull in classified peptides from the different searches
search_conditions <- list("downstream_search", "upstream_search")

for (condition in search_conditions) {
  directory <- paste0("../second_pass_", condition)
  file_path <- paste0(directory, "/analysis/output_tables/classified_peptides_data.tsv")
  loop_tbl <- read_tsv(file_path, col_types = cols("prev_aa" = col_character(),
                                                   "next_aa" = col_character(),
                                                   "frame_transition" = col_character(),
                                                   "charges" = col_character())) |> 
    rename(enz_act = source) |> 
    mutate(enz_act = str_remove(enz_act, "_pepcom"))
  assign(paste0(condition, "_classified_peptides"), loop_tbl, envir = .GlobalEnv)
  rm(directory, file_path, loop_tbl, condition)
}

all_secondpass_classified_peptides <- bind_rows(upstream = upstream_search_classified_peptides,
                                            downstream = downstream_search_classified_peptides,
                                            .id = "source") |> 
  select(!search_side)

