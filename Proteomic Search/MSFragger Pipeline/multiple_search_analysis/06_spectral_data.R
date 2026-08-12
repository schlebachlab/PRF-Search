# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

# table of all hit peptides, containing the search side and enzyme/activation method
hits_both_streams <- bind_rows(
  upstream = upstream_search_hits,
  downstream = downstream_search_hits,
  .id = "search_side"
) |> 
  select(-starts_with("n_")) |> 
  relocate(search_side, .before = enz_act)

# table that matches relevant combined_peptide.tsv files from the search side and enzyme/activation methods
hit_pepcom <- hits_both_streams |> 
  select(search_side, enz_act) |> 
  distinct() |> 
  mutate(file_path = paste0(
    "../second_pass_", search_side, "_search/results/", enz_act, "/combined_peptide.tsv"
  ))

hit_pepcom_filepaths <- hit_pepcom$file_path # character list of file paths to loop over

# make cleanup character vector for the cell line directories
cleanup_strings <- c(
  "he_la" = "HeLa",
  "gm12878" = "GM12878",
  "huvec" = "HUVEC",
  "hep_g2" = "HepG2",
  "h_esc" = "hESC",
  "k562" = "K562",
  "asp_n_hcd" = "AspN_HCD",
  "asp_n_etd" = "AspN_ETD",
  "chymo_cad" = "Chymo_CAD",
  "chymo_hcd" = "Chymo_HCD",
  "glu_c_hcd" = "GluC_HCD",
  "glu_c_etd" = "GluC_ETD",
  "lys_c_hcd" = "LysC_HCD",
  "lys_c_etd" = "LysC_ETD",
  "lys_n_etd" = "LysN_ETD",
  "lys_n_hcd" = "LysN_HCD",
  "trypsin_hcd" = "Trypsin_HCD"
) 

peptide_spectral_data <- tibble() # empty tibble to accumulate results

for (i in seq_along(hit_pepcom_filepaths)) {
  # get relevant info for each iteration
  info <- hit_pepcom |> 
    filter(file_path == hit_pepcom_filepaths[i])
  
  # temporary combined peptide tibble to get the cell lines that had spectral matches for hit peptides
  pepcom <- read_tsv(hit_pepcom_filepaths[i]) |> 
    janitor::clean_names() |> 
    mutate(search_side = info$search_side,
           enz_act = info$enz_act)|> 
    semi_join(hits_both_streams,
              join_by(search_side, enz_act, peptide_sequence == peptide))|> 
    pivot_longer(
      cols = c(ends_with("spectral_count")),
      names_to = "cell_line",
      values_to = "spectral_count"
    ) |>
    select(peptide_sequence, search_side:spectral_count) |>
    mutate(cell_line = str_remove(cell_line, "_spectral_count"),
           cell_line = str_replace_all(cell_line, cleanup_strings)) |> 
    filter(spectral_count > 0) 
  
  # append the relevant cell line info into the accumulator tibble 
  peptide_spectral_data <- bind_rows(peptide_spectral_data, pepcom)
  
  # cleanup global environment
  rm(info)
  rm(pepcom)
}

# add the relevant .tsv file paths
peptide_spectral_data <- peptide_spectral_data |> 
  mutate(
    ion_tsv = paste0(
    "../second_pass_", search_side, "_search/results/", enz_act, "/", cell_line, "/ion.tsv"),
    peptide_tsv = paste0(
      "../second_pass_", search_side, "_search/results/", enz_act, "/", cell_line, "/peptide.tsv"),
    psm_tsv = paste0(
      "../second_pass_", search_side, "_search/results/", enz_act, "/", cell_line, "/psm.tsv")
    )

# tibble of distinct relevant file paths
spectral_file_paths <- peptide_spectral_data |> 
  select(!peptide_sequence & !spectral_count) |> 
  distinct(search_side, enz_act, cell_line, .keep_all = TRUE)

# start of the for loop too pull in relevant ion.tsv files and join with peptide and protein info from database headers
hit_ion_data <- tibble() # create an empty tibble to add to for each iteration of for loop
hit_ion_tsv_files <- spectral_file_paths$ion_tsv # create a character list of ion.tsv file paths

# for loop that reads the ion.tsv, adds the cell line, search side, and enzyme activation method data
# then filters with the table of hits from both streams and joins with the same tibble to capture fasta database header info
# finally it adds the result to the empty tibble created above
for (i in seq_along(hit_ion_tsv_files)) {
  ion_tsv <- read_tsv(hit_ion_tsv_files[i]) |> 
    janitor::clean_names() |> 
    mutate(peptide = peptide_sequence,
           search_side = spectral_file_paths$search_side[i],
           enz_act = spectral_file_paths$enz_act[i],
           cell_line = spectral_file_paths$cell_line[i]) |> 
    select(!protein:mapped_proteins & !peptide_sequence) |> 
    relocate(search_side:cell_line, .before = 1) |> 
    relocate(peptide, .after = cell_line) |> 
    semi_join(hits_both_streams) |> 
    left_join(hits_both_streams) |> 
    mutate(cell_line = sapply(str_split(cell_line, "_"), function(x) paste0(x[1], collapse = "_")))
  hit_ion_data <-  bind_rows(hit_ion_data, ion_tsv)
  rm(ion_tsv)
}

# start of the for loop too pull in relevant peptide.tsv files and join with peptide and protein info from database headers
hit_peptide_data <- tibble() # create an empty tibble to add to for each iteration of for loop
hit_peptide_tsv_files <- spectral_file_paths$peptide_tsv # create a character list of peptide.tsv file paths

# for loop that reads the peptide.tsv, adds the cell line, search side, and enzyme activation method data
# then filters with the table of hits from both streams and joins with the same tibble to capture fasta database header info
# finally it adds the result to the empty tibble created above
for (i in seq_along(hit_peptide_tsv_files)) {
  peptide_tsv <- read_tsv(hit_peptide_tsv_files[i],
                          col_types = cols("Prev AA" = col_character(), 
                                           "Next AA" = col_character(), 
                                           "Charges" = col_character(),
                                           "Assigned Modifications" = col_character(),
                                           "Observed Modifications" = col_character())) |> 
    janitor::clean_names() |> 
    mutate(search_side = spectral_file_paths$search_side[i],
           enz_act = spectral_file_paths$enz_act[i],
           cell_line = spectral_file_paths$cell_line[i]) |> 
    select(!protein:mapped_proteins) |> 
    relocate(search_side:cell_line, .before = peptide) |> 
    semi_join(hits_both_streams) |> 
    left_join(hits_both_streams) |> 
    mutate(cell_line = sapply(str_split(cell_line, "_"), function(x) paste0(x[1], collapse = "_")))
  
  hit_peptide_data <- bind_rows(hit_peptide_data, peptide_tsv)
  rm(peptide_tsv)
}

# start of the for loop too pull in relevant psm.tsv files and join with peptide and protein info from database headers
hit_psm_data <- tibble() # create an empty tibble to add to for each iteration of for loop
hit_psm_tsv_files <- spectral_file_paths$psm_tsv # create a character list of ion.tsv file paths

# for loop that reads the psm.tsv, adds the cell line, search side, and enzyme activation method data
# then filters with the table of hits from both streams and joins with the same tibble to capture fasta database header info
# finally it adds the result to the empty tibble created above
for (i in seq_along(hit_psm_tsv_files)) {
  psm_tsv <- read_tsv(hit_psm_tsv_files[i]) |> 
    janitor::clean_names() |> 
    select(!protein:mapped_proteins) |> 
    relocate(peptide, .before = 1) |> 
    mutate(search_side = spectral_file_paths$search_side[i],
           enz_act = spectral_file_paths$enz_act[i],
           cell_line = spectral_file_paths$cell_line[i],
           cell_line = sapply(str_split(cell_line, "_"), function(x) paste0(x[1], collapse = "_"))) |>
    relocate(search_side:cell_line, .before = peptide) |> 
    semi_join(hits_both_streams) |> 
    left_join(hits_both_streams)
  
  hit_psm_data <- bind_rows(hit_psm_data, psm_tsv)
  rm(psm_tsv)
}

