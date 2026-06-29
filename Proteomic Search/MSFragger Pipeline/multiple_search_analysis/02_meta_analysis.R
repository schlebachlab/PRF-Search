# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

# tibble of peptides from purely native sequences found in the second pass 
# (all came from downstream and upstream search)
no_frameshift_secondpass_peptides <- all_search_classified_peptides |> 
  filter(frame_transition == "no_frameshift")

# tibble of peptides from frameshifted sequences found in the second pass 
# (found in all searches but also contains zero_frame peptides)
frameshifted_second_pass_peptides <- all_search_classified_peptides |> 
  filter(frame_transition != "no_frameshift")

# filtering join to generate tibble containing the peptides from canonical sequences
# that were also found in the first pass
native_nonframeshift_peptides_bothpasses <- no_frameshift_secondpass_peptides |> 
  rename(gene = gene_name,
         peptide_sequence = peptide) |> 
  semi_join(first_pass_peptides_and_genes, join_by(enz_act, peptide_sequence, gene))

# filtering join to generate tibble containing the peptides from canonical sequences
# that were not found in the first pass
native_nonframeshift_peptides_secondpass <- no_frameshift_secondpass_peptides |> 
  rename(gene = gene_name,
         peptide_sequence = peptide) |> 
  anti_join(first_pass_peptides_and_genes, join_by(enz_act, peptide_sequence, gene))

# summary breakdown of the second pass searches
# breakdown by search strategy, enzyme/activation, and peptide classification:
summary_total_breakdown <- all_search_classified_peptides |> 
  group_by(source, enz_act, category) |> 
  summarise(n = n())

# breakdown by search strategy, and peptide classification:
summary_search_and_pepclass <- all_search_classified_peptides |> 
  group_by(source, category) |> 
  summarise(n = n())

# breakdown by enzyme/activation, and peptide classification:
summary_enzymeactivation_and_pepclass <- all_search_classified_peptides |> 
  group_by(enz_act, category) |> 
  summarise(n = n())

# comparing what frameshift and transitional peptides are shared between
# prf_search/downstream and prf_search/upstream
main_search_hits <- main_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame")) # prf search hits

prf_upstream_shared_nonzeroframe <- main_search_hits |> 
  semi_join(upstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name)) # shared hits between prf/upstream searches (ms method, gene, peptide sequence)

test <- main_search_hits |> 
  semi_join(upstream_search_hits, join_by(enz_act, peptide, gene_name))

test2 <- upstream_search_classified_peptides |> 
  semi_join(main_search_hits, join_by(enz_act, peptide, gene_name))

test2 |> 
  group_by(category) |> 
  summarise(n = n())

prf_upstream_shared_nonzeroframe_exact <- main_search_hits |> 
  semi_join(upstream_search_hits)

upstream_unique_nonzeroframe <- upstream_search_hits |> 
  anti_join(main_search_hits, 
            join_by(enz_act, peptide, gene_name))

prf_downstream_shared_nonzeroframe <- main_search_hits |>  
  semi_join(downstream_search_classified_peptides, join_by(enz_act, peptide, gene_name)) # shared hits between prf/upstream searches

upstream_search_hits <- upstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

downstream_search_hits <- downstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

upstream_downstream_shared_nonzeroframe <- upstream_search_hits |> 
  semi_join(downstream_search_hits, join_by(enz_act, peptide, gene_name))
  



