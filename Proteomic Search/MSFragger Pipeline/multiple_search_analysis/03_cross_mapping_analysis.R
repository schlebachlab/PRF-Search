# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

# read in the cross-map files
for (condition in search_conditions) {
  directory <- paste0("../second_pass_", condition)
  upstream_path <- paste0(directory, 
                      "/analysis/possible_frameshift_peptides_classed_upstream_map.tsv")
  downstream_path <- paste0(directory, 
                            "/analysis/possible_frameshift_peptides_classed_downstream_map.tsv")
  if (file.exists(upstream_path)) {
    file_path <- upstream_path
  } else {
    file_path <- downstream_path
  }
  loop_tbl <- read_tsv(file_path)
  assign(paste0(condition, 
                if (file_path == upstream_path) {"_upstream_map"}
                else {"_downstream_map"}), 
         loop_tbl, envir = .GlobalEnv)
  rm(directory, file_path, loop_tbl, condition)
}

# classifying the main search identifications
# with the classifications from the upstream mapping file
main_search_upstream_map_classified_peptides <- main_search_classified_peptides |> 
  select(!category & !n_frameshift:n_total) |> 
  right_join(main_search_upstream_map, join_by(peptide)) |> 
  relocate(category, .after = ensembl_gene_id)

# anti-join using bottom up method, peptide sequence, gene name, and classification
# this shows the peptides that were misclassified in the initial(main) search
# with the classifications of the upstreaming mapping file
misclassified_prf_search_upstream_map <- main_search_upstream_map_classified_peptides |> 
  anti_join(main_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(upstream_category = category) |> 
  select(-starts_with("n_"))

# anti-join using bottom up method, peptide sequence, gene name, and classification
# this shows the peptides that were misclassified in the initial(main) search
# with the classifications of the initial prf mapping file
misclassified_prf_search_prf_map <- main_search_classified_peptides |> 
  anti_join(main_search_upstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(initial_prf_category = category) |> 
  select(-starts_with("n_"))

# comparison of the upstream mapping and initial prf mapping in one file
upstream_mapping_vs_prf_mapping <- misclassified_prf_search_prf_map |> 
  right_join(misclassified_prf_search_upstream_map) |> 
  relocate(upstream_category, .after = initial_prf_category)

# same as the two misclassification tibbles above
# but includes data from utility scripts about how many times a peptide fits a particular classification
misclassified_prf_search_upstream_map_with_totals <- main_search_upstream_map_classified_peptides |> 
  anti_join(main_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(upstream_category = category)

misclassified_prf_search_prf_map_with_totals <- main_search_classified_peptides |> 
  anti_join(main_search_upstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(initial_prf_category = category)

write_tsv(upstream_mapping_vs_prf_mapping, "./upstream_mapping_vs_initial_mapping.tsv")
write_tsv(misclassified_prf_search_prf_map_with_totals, "./initial_prf_with_initial_mapping.tsv")
write_tsv(misclassified_prf_search_upstream_map_with_totals, "./initial_prf_with_upstream_mapping.tsv")

# looking at the hits (transitional/peptide) from cross-mapping classification
# these are peptides classified as transitional/frameshift when:
# the initial search and downstream search results were mapped with the upstream mapping file
# and the upstream search results were mapped with the downstream mapping file

# classifying the downstream search identifications
# with the classifications from the upstream mapping file
downstream_search_upstream_map_classified_peptides <- downstream_search_classified_peptides |> 
  select(!category & !n_frameshift:n_total) |> 
  right_join(downstream_search_upstream_map, join_by(peptide)) |> 
  relocate(category, .after = ensembl_gene_id)

# classifying the upstream search identifications
# with the classifications from the downstream mapping file
upstream_search_downstream_map_classified_peptides <- upstream_search_classified_peptides |> 
  select(!category & !n_frameshift:n_total) |> 
  right_join(upstream_search_downstream_map, join_by(peptide)) |> 
  relocate(category, .after = ensembl_gene_id)

# filter the cross-mapped peptides for hits (transition/frameshift peptides)
downstream_search_upstream_map_hits <- downstream_search_upstream_map_classified_peptides |> 
  filter(!str_detect(category, "zero_frame") & !str_detect(category, "not_found"))

upstream_search_downstream_map_hits <- upstream_search_downstream_map_classified_peptides |> 
  filter(!str_detect(category, "zero_frame") & !str_detect(category, "not_found"))

main_search_upstream_map_hits <- main_search_upstream_map_classified_peptides |> 
  filter(!str_detect(category, "zero_frame") & !str_detect(category, "not_found"))

# filter upstream / downstream cross-mappings for not_found peptides
downstream_search_upstream_map_not_found <- downstream_search_upstream_map_classified_peptides |> 
  filter(str_detect(category, "not_found"))

upstream_search_downstream_map_not_found <- upstream_search_downstream_map_classified_peptides |> 
  filter(str_detect(category, "not_found"))

