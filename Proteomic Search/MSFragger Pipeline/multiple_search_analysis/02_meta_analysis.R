# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

library(tidyverse)

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
# prf_search/downstream and prf_search/upstream:

# filter out peptides with any zero_frame classification from each search
main_search_hits <- main_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame")) # prf search hits

downstream_search_hits <- downstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

upstream_search_hits <- upstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

# begin comparison of each search:

# shared hits between prf/upstream searches (ms method, gene, peptide sequence)
# NOTE: because the main search hits were filtered with the upstream search hits,
# classifications are based on the mapping file of the original search
prf_upstream_shared_nonzeroframe <- main_search_hits |> 
  semi_join(upstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# shared hits between prf/upstream searches (exact matches)
# the difference between this and above is different classification based on different mapping files
# NOTE: because the main search hits were filtered with the upstream search hits,
# classifications are based on the mapping file of the original search
prf_upstream_shared_nonzeroframe_exact <- main_search_hits |> 
  semi_join(upstream_search_hits)

# unique transitional/frameshift peptides found in the upstream search
# NOTE: because the upstream search hits were filtered by the main search hits,
# classifications are based on the mapping file of the main search
upstream_unique_nonzeroframe <- upstream_search_hits |> 
  anti_join(main_search_hits, 
            join_by(enz_act, peptide, gene_name))

# shared identification between prf/downstream searches (ms method, gene, peptide sequence)
# NOTE: because the main search hits were filtered with the downstream search identifications (not hits),
# classifications are based on the mapping file of the main search
# NOTE: doing a filter with exact matches (classifications) produces zero results,
# this is because the classification of the downstream search with its respective mapping file
# will be very different
prf_downstream_ids_shared_nonzeroframe <- main_search_hits |>  
  semi_join(downstream_search_classified_peptides, join_by(enz_act, peptide, gene_name))
# this shows the hit peptides from the main search that match to identified peptides in the downstream search
# returns 2 peptides

prf_downstream_hits_shared_nonzeroframe <- main_search_hits |>  
  semi_join(downstream_search_hits, join_by(enz_act, peptide, gene_name))
# this shows the hit peptides from the main search that were also hits in the downstream search
# returns 1 peptide - the missing peptide from above was classified as mixture_with_zero_frame,
# when classification was done with the downstream mapping file

# shared hits between upstream/downstream searches (ms method, gene, peptide sequence)
# NOTE: because the upstream search hits were filtered with the downstream search hits,
# classifications are based on the mapping file of the upstream search
upstream_downstream_shared_nonzeroframe <- upstream_search_hits |> 
  semi_join(downstream_search_hits, join_by(enz_act, peptide, gene_name))

prf_search_hits_in_upstream_ids <- main_search_hits |> 
  semi_join(upstream_search_hits, join_by(enz_act, peptide, gene_name)) # this returns 12 hits
# this is the main search hits filtered using the hits from the upstream search
# 12 hits classified as transitional/frameshift in the upstream search, 
# when mapped to the upstream search mapping file
# 6 hits in the upstream search are unique to that search

upstream_identifications_in_prf_hits <- upstream_search_classified_peptides |> 
  semi_join(main_search_hits, join_by(enz_act, peptide, gene_name)) # this returns 25 hits
# this is all of the upstream search identifications (mapped with the upstream mapping file)
# these identifications are filtered with the main search hits
# there are 13 peptides that are classified as mixture_with_zero_frame, 
# these 13 peptides are classified as frameshift/transitional in the main search
# this shows that the mapping of the main search using the main search mapping file is incomplete


upstream_identifications_in_prf_hits |> 
  group_by(category) |> 
  summarise(n = n())

## keep the global environment for the cross-mapping analysis script

