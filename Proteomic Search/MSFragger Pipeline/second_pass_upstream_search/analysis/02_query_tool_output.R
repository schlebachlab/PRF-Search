# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

# after querying the peptides with fs_pep_utilities.py, use this script to analyse the output
# the global environment of the previous script is needed for this script to function

#read output into tibbles
all_peptides_query <- read_tsv("all_peptides_classed.tsv")
frameshift_peptides_query <- read_tsv("possible_frameshift_peptides_classed.tsv")

#join the second pass output with the query output to capture the meta data with the peptide classification
all_peptides_class <- all_peptides |> 
  select(source:protein) |> 
  rename(peptide = peptide_sequence) |> 
  left_join(all_peptides_query, join_by(peptide))

frameshift_peptides_class <- predicted_frameshift_peptides |> 
  select(source:protein) |> 
  rename(peptide = peptide_sequence) |> 
  left_join(frameshift_peptides_query, join_by(peptide))

# widen the data to capture the info in Chuck's headers
prdct_frmshft_pep_class_parsed <- frameshift_peptides_class |> 
  separate_rows(protein, sep = "\\|")|> 
  separate(protein, into = c("key", "value"), sep = "=")|> 
  pivot_wider(names_from = key, values_from = value) |> 
  relocate(gene_name, .after = peptide) |> 
  relocate(ensembl_transcript_id:ensembl_gene_id, .after = gene_name) |> 
  relocate(category, .after = ensembl_gene_id)

#make a tibble with the canonical proteins and contaminants
canon_contam_meta <- all_peptides_class |> 
  filter(str_detect(protein, "sp\\|"))

# write out data to .tsv files - if the files are not needed, just make the tibbles
write.table(prdct_frmshft_pep_class_parsed, file = './output_tables/classified_peptides_data.tsv', 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE, 
            sep = "\t")

transitional_peptides <- prdct_frmshft_pep_class_parsed |> 
  filter(category == 'transitional_only')

write.table(transitional_peptides, file = './output_tables/transitional_seq_peptides.tsv', 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE, 
            sep = "\t")

frameshift_only_peptides <- prdct_frmshft_pep_class_parsed |> 
  filter(category == 'frameshift_only')

write.table(frameshift_only_peptides, file = './output_tables/frameshift_seq_peptides.tsv', 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE, 
            sep = "\t")

hit_peptides <- prdct_frmshft_pep_class_parsed |> 
  filter(category == 'transitional_only' | category == 'frameshift_only')

write.table(hit_peptides, file = './output_tables/transitional_and_frameshift_peptides.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

zeroframe_frameshift_stop <- prdct_frmshft_pep_class_parsed |> 
  filter(category == "zero_frame_with_frameshift_stop_only")

write.table(zeroframe_frameshift_stop, file = './output_tables/zeroframe_frameshift_stop.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
