# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

# the global environment from previous scripts is needed for proper function of this script.

# a for loop to join the metadata in frameshifted fasta headers 
# with the MS data from the combined_peptide.tsv, with canonical uniprot IDs (ie sp|XX_YY) removed
# this returns tibbles with the parsed frameshifted fasta headers and associated MS search data
# for all peptides predicted from frameshift prediction algorithm
for (df in peptide_list_names) {
  og_tbl <- get(df)
  meta_filt <- prdct_frmshft_pep_class_parsed |> 
    filter(source == df)
  mod_tbl <- og_tbl |> 
    select(!protein:mapped_proteins) |> 
    rename(peptide = peptide_sequence) |> 
    inner_join(meta_filt)
  new_name <- paste0(df, "_meta")
  assign(new_name, mod_tbl, envir = .GlobalEnv)
  rm(og_tbl)
  rm(meta_filt)
  rm(mod_tbl)
  rm(new_name)
}

# a for loop to join the metadata in frameshifted fasta headers 
# with the MS data for the frameshifted and transitional peptides from the combined_peptide.tsv
# this returns tibbles with the parsed frameshifted fasta headers and associated MS search data
# for transitional and frameshifted peptides
for (df in peptide_list_names) {
  og_tbl <- get(df)
  hit_filt <- hit_peptides |> 
    filter(source == df)
  mod_tbl <- og_tbl |> 
    select(!protein:mapped_proteins) |> 
    rename(peptide = peptide_sequence) |> 
    inner_join(hit_filt)
  new_name <- paste0(df, "_hits")
  assign(new_name, mod_tbl, envir = .GlobalEnv)
  rm(og_tbl)
  rm(hit_filt)
  rm(mod_tbl)
  rm(new_name)
}

# a for loop to pivot the hits from each enzyme and activation method
# creates a longer tibble to identify the BU experiment where hits were ID'd
# resulting ENZ_ACT_hits is a pivot table of the transitional and frameshift peptides
# with corresponding MS experiment
for (df in peptide_list_names) {
  tbl_name <- paste0(df, "_hits")
  og_tbl <- get(tbl_name)
  tbl_long <- og_tbl |> 
    pivot_longer(
      cols = c(ends_with("spectral_count"), ends_with("intensity")),
      names_to = "ms_data_type",
      values_to = "ms_data_value"
    ) |> 
    relocate(ms_data_type:ms_data_value, .after = peptide)
  tbl_name <- str_remove(tbl_name, "_pepcom")
  assign(tbl_name, tbl_long, envir = .GlobalEnv)
  rm(og_tbl)
  rm(tbl_long)
  rm(tbl_name)
}

# create a large tibble for hits across all enzyme and activation conditions
# by combing the long pivots of ENZ_ACT_hits, filter for MS experiments with identified scans,
# and provide cell-line and enz_act info to help match to MS spectrum file
hit_list_names <- paste0(enz_act_list, "_hits") #make a list of the names of hit tibbles
cleanup_strings <- c(
  "hela" = "HeLa",
  "gm12878" = "GM12878",
  "huvec" = "HUVEC",
  "hepg2" = "HepG2",
  "hesc" = "hESC",
  "k562" = "K562",
  "aspn_hcd" = "AspN_HCD",
  "aspn_etd" = "AspN_ETD",
  "chymo_cad" = "Chymo_CAD",
  "chymo_hcd" = "Chymo_HCD",
  "gluc_hcd" = "GluC_HCD",
  "gluc_etd" = "GluC_ETD",
  "lysc_hcd" = "LysC_HCD",
  "lysc_etd" = "LysC_ETD",
  "lysn_etd" = "LysN_ETD",
  "lysn_hcd" = "LysN_HCD",
  "trypsin_hcd" = "Trypsin_HCD"
) #make cleanup character vector

#make a table of all identified hits with the cell line, enzyme, and activation method decoded to set up parsing psm.tsv
hit_list_notryp <- mget(hit_list_names) |> 
  bind_rows(.id = "source") |> 
  filter(source != "Trypsin_CAD_hits") |> 
  select(!ends_with("match_type")) |> 
  filter(ms_data_value > 0) |> 
  relocate(source, .before = 1) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "he_la", "hela")) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "h_esc", "hesc")) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "hep_g2", "hepg2")) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "lys_c", "lysc")) |>
  mutate(ms_data_type = str_replace(ms_data_type, "lys_n", "lysn")) |>
  mutate(ms_data_type = str_replace(ms_data_type, "glu_c", "gluc")) |>
  mutate(ms_data_type = str_replace(ms_data_type, "asp_n", "aspn")) |>
  mutate(
    cell_line = sapply(str_split(ms_data_type, "_"), function(x) paste0(x[1], collapse = "_")),
    .before = ms_data_type
  ) |> 
  mutate(
    enz_act = sapply(str_split(ms_data_type, "_"), function(x) paste0(x[2:3], collapse = "_")),
    .before = ms_data_type
  ) |> 
  mutate(source = str_remove(source, "_hits")) |> 
  rename(exp_group = source) |> 
  mutate(cell_line = str_replace_all(cell_line, cleanup_strings)) |> 
  mutate(enz_act = str_replace_all(enz_act, cleanup_strings)) |> 
  relocate(gene_name, .before = 1)

#since the file structure of the Tryp_CAD second pass search is slightly different 
#(Due to the manifest experiment group names) a separate hit list needs to be prepared
Trypsin_CAD_hit_list <- Trypsin_CAD_hits |> 
  select(!ends_with("match_type")) |> 
  filter(ms_data_value > 0) |> 
  relocate(source, .before = 1) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "he_la", "hela")) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "h_esc", "hesc")) |> 
  mutate(ms_data_type = str_replace(ms_data_type, "hep_g2", "hepg2")) |> 
  mutate(
    cell_line = sapply(str_split(ms_data_type, "_"), function(x) paste0(x[1], collapse = "_")),
    .before = ms_data_type
  ) |> 
  mutate(enz_act = "Trypsin_HCD", .after = cell_line) |> 
  mutate(source = str_remove(source, "_pepcom")) |> 
  rename(exp_group = source) |> 
  mutate(cell_line = str_replace_all(cell_line, cleanup_strings))

hit_list <- hit_list_notryp |> 
  bind_rows(Trypsin_CAD_hit_list)

# make a tibble with the frameshift and transitional peptides with a column specifying the file path
# for the relevant psm.tsv
#first for Trypsin:
Trypsin_cell_exp <- Trypsin_CAD_hit_list |> 
  select(exp_group:enz_act) |> 
  mutate(file_path = paste0("../results/", exp_group, "/", cell_line, "/psm.tsv")) |> 
  distinct(exp_group, peptide, file_path, .keep_all = TRUE)
#then for the rest, binding to tryspin at the end:
pep_cell_exp <- hit_list_notryp |> 
  select(exp_group:enz_act) |> 
  mutate(file_path = paste0("../results/", exp_group, "/" ,cell_line, "_", enz_act, "/psm.tsv")) |> 
  distinct(exp_group, peptide, file_path, .keep_all = TRUE) |> 
  bind_rows(Trypsin_cell_exp)

# iterate over the list of peptides and pull the relvant psm.tsv files into tibbles
psm_tibble_list <- list() #create an empty list for the following for loop to paste into
#a for loop that iterates over the rows in the pep_cell_exp tibble, grabbing the experiment group, cell line, and activation
#it pulls the psm.tsv and then names it based on those parameters, then adds it to the tibble list for subsequent code

relevant_fragpipe_files <- pep_cell_exp |> 
  select(!peptide) |> 
  distinct(exp_group, cell_line, enz_act, file_path)

for (i in seq_len(nrow(relevant_fragpipe_files))) {
  exp <- relevant_fragpipe_files$exp_group[i]
  cell <- relevant_fragpipe_files$cell_line[i]
  enz_act <- relevant_fragpipe_files$enz_act[i]
  file_path <- relevant_fragpipe_files$file_path[i]
  df <- read_tsv(file_path)
  assign(paste0(exp, "_", cell, "_", enz_act), df) |> 
    janitor::clean_names()
  psm_tibble_list[[i]] <- paste0(exp, "_", cell, "_", enz_act)
  rm(df)
} 

# pool the psm data into one tibble and clean up
psm_tibble_list <- unique(psm_tibble_list) |> 
  unlist()
relevant_psm <- mget(psm_tibble_list, envir = .GlobalEnv) |> 
  bind_rows(.id = "source") |> 
  janitor::clean_names() |> 
  mutate(exp_group = sapply(str_split(source, "_"), function(x) paste0(x[1:2], collapse = "_")), 
    .after = 1) |> 
  mutate(cell_line = sapply(str_split(source, "_"), function(x) paste0(x[3])),
    .after = 2) |>
  mutate(enz_act = sapply(str_split(source, "_"), function(x) paste0(x[4:5], collapse = "_")),
    .after = 3) |> 
  select(!protein:mapped_proteins)

# filter the hit list for only spectral count and then left join to the psm tibble with scan info
hit_list_speccount <- hit_list |> 
  filter(str_detect(ms_data_type, "spectral_count"))

hit_peptides_spectral_data <- hit_list_speccount |> 
  left_join(relevant_psm) |> 
  relocate(spectrum:spectrum_file, .after = ms_data_value) |> 
  relocate(ensembl_transcript_id:category, .after = gene_name) |> 
  relocate(ensembl_uniprot_annotation:blastp_uniprot_annotation, .after = category) |> 
  relocate(ensembl_transcript_id:blastp_uniprot_annotation, .after = spectrum_file)

# write tables
write.table(hit_list, file = './output_tables/transitional_and_frameshift_peptides_MSdata.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(hit_peptides_spectral_data, file = './output_tables/hit_peptides_spectral_data.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
