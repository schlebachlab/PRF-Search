Gene Ontology Analysis

This folder contains the code used to perform Gene Ontology (GO) enrichment analyses on subsets of TMD-slip (Harrington) motifs. The script consumes a previously generated Harrington motif table (from the main motif-search pipeline), maps motif-associated transcripts to UniProt accessions, and tests for GO term enrichment using a UniProt-derived GO annotation background.

The analysis focuses on reproducing the manuscript’s GO enrichment results while keeping motif definitions consistent with the discovery pipeline. The theoretical rationale for motif definitions and subset choices is described in the manuscript and Supplemental Theory; this README focuses on running the script and understanding its outputs.

The primary script in this folder is gene_ontology_enrichment_scheme_goat__no_stemloop_required_35-55_gap.py, which is run in configuration-file mode.

Running the GO enrichment pipeline

Run the script with the configuration file:

python gene_ontology_enrichment_scheme_goat__no_stemloop_required_35-55_gap.py config_file_v84_gene_ontology.ini

All behavior is controlled by the configuration file; no additional command-line arguments are required.

Configuration file

The manuscript analysis uses a config file with these fields:

[Paths]
path_motifs_csv = motif_search_output_v73_motif_search_harrington_motifs.csv
path_uniprot_tsv = uniprotkb_AND_reviewed_true_AND_model_o_2024_10_02.tsv
blastp_results_path_csv = motif_search_output_v73_motif_search_ensembl_uniprot_dict_full_blastp.csv
path_tmd_csv = topcons_ensembl_output_march_2025.csv
gene_ontology_obo_path = go-basic.obo

[Parameters]
outname = v84_no_stemloop_required_gap_35-55

The key input is path_motifs_csv, which should be the Harrington motif table produced by the motif-search pipeline. The UniProt TSV is expected to contain reviewed human proteins with GO term annotations. The GO ontology file should be the standard go-basic.obo.

What the script does

The script loads the Harrington motif table and splits motifs into commonly used motif subsets (strict vs loose slip sites; canonical vs alternative isoforms; and “uniquely alternative” motifs). For each subset it collects UniProt accessions, defines an appropriate UniProt background set, and performs GO enrichment analyses across the three GO namespaces (Biological Process, Cellular Component, and Molecular Function).

In addition to GO terms, the script performs keyword enrichment testing for the UniProt keyword “Disease variant” across the same motif subsets.

Outputs

All output files are prefixed using the value of outname.

The main output is a multi-sheet Excel workbook:

outname_gene_ontology_results.xlsx

This workbook contains one sheet per motif subset and GO namespace, and includes enrichment statistics (odds ratios, log-odds, fold enrichment, raw p-values, and FDR-adjusted p-values) along with counts used to form the enrichment contingency tables.

The script also emits additional supporting outputs used for bookkeeping and downstream inspection:

outname_motifs_spacing_only_disease_and_gpcr_log.txt
A log file summarizing counts of motifs and motif-associated genes with UniProt keywords such as “Disease variant,” “G-protein coupled receptor,” and “Ion channel,” stratified across motif subsets.

outname_motifs_unique_to_alternative_isoforms.csv
A CSV table listing motifs that appear uniquely in alternative isoforms (under the subset definition used in the script).

outname_motif_table_with_new_annotations.csv
A version of the motif table with additional annotations appended (e.g., GPCR / ion channel / disease keyword flags and UniProt keyword strings), intended for manual inspection and figure generation.
