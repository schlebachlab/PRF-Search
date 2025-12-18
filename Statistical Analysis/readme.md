Statistical analysis of TMD-slip motifs

This folder contains the scripts used to perform statistical analyses on TMD-slip motifs identified in the transcriptome. These analyses quantify enrichment, positional bias, isoform usage, and spacing constraints, and they generate the statistics reported in the manuscript figures and tables.

The statistical analysis pipeline reuses much of the transcript parsing, TMD association, and slip-site scanning logic from the motif search pipeline, but it differs in two key ways. First, no frameshift translation products or chimeric FASTA files are generated. Second, the script performs a series of explicit statistical tests on motif distributions rather than emitting peptide-level outputs.

The primary script in this folder is enhanced_harrington_search_v84_one_hept_stats_first_set.py, which is run in configuration-file mode for reproducibility.

Running the statistical analysis

Run the script with the configuration file:

python enhanced_harrington_search_v84_one_hept_stats_first_set.py config_file_v84_one_hept_stats_first_set.ini

All analysis behavior is controlled by the configuration file; no additional command-line arguments are required.

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
outname = motif_search_output_v84_one_hept_stats_first_set
scoring_sys = 2

These inputs mirror those used in the motif search pipeline, ensuring that motif definitions and background transcript sets are consistent across analyses.

Overview of analyses performed

The script identifies candidate TMD-slip motifs using the same spacing and friction criteria as the motif search pipeline, but instead of emitting translation products it aggregates motif-level and transcript-level statistics.

Analyses performed include:

• Enumeration of canonical versus alternative transcripts in the Ensembl CDS database and among TMD-containing transcripts
• Counts of transcripts, genes, and HGNC symbols with predicted transmembrane domains
• Distribution of slip-site distances relative to upstream TMDs
• Comparison of strict versus loose slippery heptamer classes
• Enrichment of motifs at the ideal TMD-to-slip-site spacing
• Background expectations derived from all candidate heptamers in TMD-containing transcripts

Several statistical tests are applied, including nonparametric distribution comparisons, contingency-table analyses, and binomial tests. These tests are implemented directly in the script and use SciPy for statistical evaluation.

Outputs

All output files are prefixed using the value of outname.

Primary outputs include:

• A CSV file containing all identified TMD-slip motifs with full transcript- and TMD-level annotations
• CSV files summarizing slip-site–to–TMD distance distributions for strict and loose slippery heptamers
• Log files reporting transcript counts, gene counts, isoform breakdowns, and motif tallies
• Summary tables used to generate the statistical figures reported in the manuscript

Because no frameshift products are generated in this pipeline, no FASTA files or peptide-level mapping tables are produced.
