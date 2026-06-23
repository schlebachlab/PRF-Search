# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

# list of enzyme and activation method
enz_act_list <- c("AspN_CAD", "AspN_ETD", "Chymo_CAD", "Chymo_ETD", "GluC_CAD", "GluC_ETD", "LysC_CAD", "LysC_ETD", "LysN_CAD", "LysN_ETD", "Trypsin_CAD")

# loop over each method and pull .tsv into a tibble
for (method in enz_act_list) {
  file_path <- paste0("../results/", method, "/combined_peptide.tsv")
  pepcom <- read_tsv(file_path, col_types = cols("Prev AA" = col_character(), "Charges" = col_character()))
  pepcom <- janitor::clean_names(pepcom)
  pepcom$prev_aa <- as.character(pepcom$prev_aa)
  pepcom$next_aa <- as.character(pepcom$next_aa)
  assign(paste0(method, "_pepcom"), pepcom)
  rm(pepcom)
}

# make a list of the tibbles created in the for loop and combine them
peptide_list_names <- paste0(enz_act_list, "_pepcom")
peptide_lists <- mget(peptide_list_names, envir = .GlobalEnv)
all_peptides <- bind_rows(peptide_lists, .id = "source") |> 
  select(source:mapped_proteins)

# filter peptide list for frameshifted peptides and write to .tsv and .txt file
predicted_frameshift_peptides <- all_peptides |> 
  filter(str_detect(protein, "ensembl_gene_id"))

frmshift_pep <- predicted_frameshift_peptides |> 
  select(source:protein)
write.table(predicted_frameshift_peptides, file = 'possible_frameshift_peptides.tsv', row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
frmshift_pep_list <- predicted_frameshift_peptides |> 
  select(peptide_sequence)
write.table(frmshift_pep_list, file = 'possible_frameshift_peptides.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# write .txt and .tsv files of all identified peptides
peptide_ids <- all_peptides |> 
  select(source:protein)
write.table(peptide_ids, file = 'all_peptides.tsv', row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
peptide_ids_list <- peptide_ids |> 
  select(peptide_sequence)
write.table(peptide_ids_list, file = 'all_peptides.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

#end of this script; keep the global environment from this script and use Chuck's query tool on outputs
