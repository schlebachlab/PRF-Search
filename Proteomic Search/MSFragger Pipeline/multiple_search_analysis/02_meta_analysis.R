# Author: Ali Farzam (GitHub: theAliFarzam)
# Correspondence: afarzam@purdue.edu | Dr. Bryon S. Drown (bsdrown@purdue.edu)

# summary breakdown of the second pass searches
# breakdown by search strategy, enzyme/activation, and peptide classification:
summary_total_breakdown <- all_secondpass_classified_peptides |> 
  group_by(source, enz_act, category) |> 
  summarise(n = n())

# breakdown by search strategy, and peptide classification:
summary_search_and_pepclass <- all_secondpass_classified_peptides |> 
  group_by(source, category) |> 
  summarise(n = n())

# breakdown by enzyme/activation, and peptide classification:
summary_enzymeactivation_and_pepclass <- all_secondpass_classified_peptides |> 
  group_by(enz_act, category) |> 
  summarise(n = n())

# breakdown of peptide categorization for identifications made in both passes
summary_peptide_category_both_passes <- peptides_bothpasses |> 
  group_by(category) |> 
  summarise(n = n())

# breakdown of peptide categorization for identifications made in the second pass
summary_peptide_category_secondpass <- peptides_secondpass_only |> 
  group_by(category) |> 
  summarise(n = n())

# comparing what frameshift and transitional peptides are shared between
# downstream and /upstream:

# filter out peptides with any zero_frame classification from each search
downstream_search_hits <- downstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

upstream_search_hits <- upstream_search_classified_peptides |> 
  filter(!str_detect(category, "zero_frame"))

# begin comparison of each search:

# transitional and frameshift peptides identified in the upstream search
# also found in the downstream search
upstream_hits_in_downstream_search <- upstream_search_hits |> 
  semi_join(downstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# peptide identifications found in the upstream search
# also found in the downstream search
upstream_ids_in_downstream_search <- upstream_search_classified_peptides |> 
  semi_join(downstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# unique transitional/frameshift peptides found in the upstream search
upstream_unique_hits <- upstream_search_hits |> 
  anti_join(downstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# transitional and frameshift peptides identified in the downstream search
# also found in the upstream search
downstream_hits_in_upstream_search <- downstream_search_hits |> 
  semi_join(upstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# peptide identifications found in the downstream search
# also found in the upstream search
downstream_ids_in_upstream_search <- downstream_search_classified_peptides |> 
  semi_join(upstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

# unique transitional/frameshift peptides found in the downstream search
downstream_unique_hits <- downstream_search_hits |> 
  anti_join(upstream_search_classified_peptides, 
            join_by(enz_act, peptide, gene_name))

## keep the global environment for the cross-mapping analysis script

