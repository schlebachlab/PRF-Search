# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

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

# classifying the downstream search identifications
# with the classifications from the upstream mapping file
downstream_search_upstream_map_classified_peptides <- downstream_search_classified_peptides |> 
  select(!category & !n_frameshift:n_total) |> 
  right_join(downstream_search_upstream_map, join_by(peptide)) |> 
  relocate(category, .after = ensembl_gene_id)

# anti-join using bottom up method, peptide sequence, gene name, and classification
# this shows the peptides that were misclassified in the downstream search
# with the classifications of the upstreaming mapping file
misclassified_downstream_search_upstream_map <- downstream_search_upstream_map_classified_peptides |> 
  anti_join(downstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(upstream_category = category) |> 
  select(-starts_with("n_"))

# anti-join using bottom up method, peptide sequence, gene name, and classification
# this shows the peptides that were misclassified in the initial(main) search
# with the classifications of the initial prf mapping file
misclassified_downstream_search_downstream_map <- downstream_search_classified_peptides |> 
  anti_join(downstream_search_upstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(downstream_category = category) |> 
  select(-starts_with("n_"))

# comparison of the upstream mapping and initial prf mapping in one file
downstream_search_mapping_comparison <- misclassified_downstream_search_downstream_map |> 
  right_join(misclassified_downstream_search_upstream_map) |> 
  relocate(upstream_category, .after = downstream_category)

# same as the two misclassification tibbles above
# but includes data from utility scripts about how many times a peptide fits a particular classification
misclassified_downstream_search_upstream_map_with_totals <- downstream_search_upstream_map_classified_peptides |> 
  anti_join(downstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(upstream_category = category)

misclassified_downstream_search_downstream_map_with_totals <- downstream_search_classified_peptides |> 
  anti_join(downstream_search_upstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(downstream_category = category)

# looking at the hits (transitional/peptide) from cross-mapping classification
# these are peptides classified as transitional/frameshift when:
# the initial search and downstream search results were mapped with the upstream mapping file
# and the upstream search results were mapped with the downstream mapping file

# classifying the upstream search identifications
# with the classifications from the downstream mapping file
upstream_search_downstream_map_classified_peptides <- upstream_search_classified_peptides |> 
  select(!category & !n_frameshift:n_total) |> 
  right_join(upstream_search_downstream_map, join_by(peptide)) |> 
  relocate(category, .after = ensembl_gene_id)

reclassified_upstream_search_downstream_map <- upstream_search_downstream_map_classified_peptides |> 
  anti_join(upstream_search_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(downstream_category = category) |> 
  select(-starts_with("n_"))

reclassified_upstream_search_upstream_map <- upstream_search_classified_peptides |> 
  anti_join(upstream_search_downstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, category)) |> 
  rename(upstream_category = category) |> 
  select(-starts_with("n_"))

upstream_search_mapping_comparison <- reclassified_upstream_search_downstream_map |> 
  right_join(reclassified_upstream_search_upstream_map) |> 
  relocate(upstream_category, .after = downstream_category)

# filter the cross-mapped peptides for hits (transition/frameshift peptides)
downstream_search_upstream_map_hits <- downstream_search_upstream_map_classified_peptides |> 
  filter(!str_detect(category, "zero_frame") & !str_detect(category, "not_found"))

upstream_search_downstream_map_hits <- upstream_search_downstream_map_classified_peptides |> 
  filter(!str_detect(category, "zero_frame") & !str_detect(category, "not_found"))

# filter upstream / downstream cross-mappings for not_found peptides
downstream_search_upstream_map_not_found <- downstream_search_upstream_map_classified_peptides |> 
  filter(str_detect(category, "not_found"))

upstream_search_downstream_map_not_found <- upstream_search_downstream_map_classified_peptides |> 
  filter(str_detect(category, "not_found"))

# cross-map analysis of the hits in downstream search reclassified with upstream map
# and vice versa
downstream_hits_reclassified <- downstream_search_hits |> 
  rename(downstream_map_category = category) |> 
  select(-starts_with("n_")) |> 
  left_join(downstream_search_upstream_map_classified_peptides) |> 
  rename(upstream_map_category = category) |> 
  relocate(upstream_map_category, .after = downstream_map_category) |> 
  select(-starts_with("n_"))

downstream_hits_not_identified_in_upstream_search <- downstream_hits_reclassified |> 
  filter(upstream_map_category != "not_found")

# this is just to see how the upstream hits are reclassified with the downstream map
# the peptides that are not shared between the downstream and upstream hits are not found
upstream_hits_reclassified <- upstream_search_hits |> 
  rename(upstream_map_category = category) |> 
  select(-starts_with("n_")) |> 
  left_join(upstream_search_downstream_map_classified_peptides) |> 
  rename(downstream_map_category = category) |> 
  relocate(downstream_map_category, .before = upstream_map_category) |> 
  select(-starts_with("n_"))

# create a large tibble with relevant cross-map information:
# make a meta table of all peptide identifications with their different classifications
cross_map_meta_table <- upstream_search_classified_peptides |> 
  rename(upstream_category = category) |> 
  select(!prev_aa:terminal_stop_codon) |> 
  right_join(upstream_search_downstream_map_classified_peptides) |> 
  rename(upstream_crossmap_category = category) |> 
  relocate(upstream_crossmap_category, .after = upstream_category) |> 
  select(!prev_aa:n_total) |>  
  full_join(downstream_search_classified_peptides,
            join_by(enz_act, peptide, gene_name)) |> 
  rename(downstream_category = category) |> 
  relocate(downstream_category, .after = upstream_crossmap_category) |> 
  select(!prev_aa:terminal_stop_codon) |> 
  left_join(downstream_search_upstream_map_classified_peptides,
            join_by(enz_act, peptide, gene_name, 
                    ensembl_transcript_id.y == ensembl_transcript_id, 
                    ensembl_gene_id.y == ensembl_gene_id)) |> 
  rename(downstream_crossmap_category = category,
        ensembl_gene_id_downstream = ensembl_gene_id.y,
        ensembl_transcript_id_downstream = ensembl_transcript_id.y,
        ensembl_gene_id_upstream = ensembl_gene_id.x,
        ensembl_transcript_id_upstream = ensembl_transcript_id.x) |> 
  select(!prev_aa:n_total) |> 
  relocate(downstream_crossmap_category, .after = downstream_category) |> 
  relocate(upstream_category:downstream_crossmap_category, .after = gene_name)

# table of peptides identified in both upstream and downstream searches
# assigned to different isoforms in upstream vs downstream
fragpipe_mismatched_peptides <-  cross_map_meta_table |> 
  filter(ensembl_transcript_id_upstream != ensembl_transcript_id_downstream)



