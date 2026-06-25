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
rm(peptide_lists, envir = .GlobalEnv) #free up RAM

peptide_and_gene <- all_peptides |> 
  select(source:charges, gene)

pep_gene_chunks <- split(peptide_and_gene, ceiling(seq_len(nrow(peptide_and_gene)) / 8e5))
imap(pep_gene_chunks, ~ write_tsv(.x, paste0("first_pass_peptide_and_genes_chunk", .y, ".tsv")))
