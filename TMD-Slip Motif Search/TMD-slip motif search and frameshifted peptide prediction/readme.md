Motif search and frameshifted peptide prediction

This folder contains the main pipeline script used to identify TMD-slip motifs and to generate predicted −1 frameshift translation products and peptide libraries used for proteomics validation.

The primary script is enhanced_harrington_search_single_ss_per_tmd_v120_fix_motif_table_all.py. The pipeline is run in configuration-file mode for reproducibility.

Running the pipeline

Run the script with the configuration file:

python enhanced_harrington_search_single_ss_per_tmd_v120_fix_motif_table_all.py config_file_v120_motif_search_best_ss.ini

Configuration file

The manuscript analysis uses the following configuration:

[Paths]
path_fasta = Homo_sapiens.GRCh38.cds.all_release-113_2024-08-15.fa
path_uniprot = uniprotkb_taxonomy_id_9606_AND_reviewed_2024_09_07.xlsx
blastp_results_path = blastp_results_uniprot_all_hs_rev_ensembl_r113_all.tsv
path_tmd_csv = topcons_ensembl_output_march_2025.csv
path_slipsites = noncanonical_slipsites_2.8_cutoff.txt
path_friction_csv = all_heptamers_friction_scan.csv
path_ensembl_features = ensembl_features_12_March_2025.tsv

[Parameters]
gap_near = 35
gap_far = 55
outname = motif_search_singleSS_v120_fs_pep
scoring_sys = 2
all_ss_candidates = 2

Key outputs

All outputs are prefixed using the value of outname.

The primary motif table used in the manuscript is written as:
<outname>_harrington_motifs.csv

A broader table of all motif candidates is also written as:
<outname>_all_motifs.csv

Predicted full-length chimeric translation products are written as gzipped FASTA files:
<outname>_all_full_length_chimeras.fa.gz
<outname>_canonical_full_length_chimeras.fa.gz

The pipeline also digests chimeric products and emits peptide FASTAs and peptide-to-parent mapping tables used in downstream proteomics analyses. These include:
<outname>_digested_fs_pep_all.fa.gz
<outname>_peptide2parent_all.tsv.gz
<outname>_digested_fs_pep_canonical.fa.gz
<outname>_peptide2parent_canonical.tsv.gz

Because these peptide-level outputs can be very large, they are typically regenerated locally and not tracked directly in version control.
