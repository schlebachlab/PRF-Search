Transmembrane domain analysis

This folder contains utilities used to generate refined transmembrane-domain (TMD) annotations from Ensembl CDS transcripts. The outputs of this stage are required inputs for the TMD-slip motif search.

Translating Ensembl CDS transcripts

The script translate_cds_v2.py translates Ensembl CDS nucleotide sequences into protein sequences while preserving Ensembl transcript identifiers. This is used instead of Ensembl-provided peptide FASTAs to ensure consistent transcript-level identifiers throughout the pipeline.

Splitting FASTA files for TOPCONS2 submission

Because the full translated FASTA may be too large for submission to the TOPCONS2 web server, the utility split_fasta.py can be used to split the FASTA into smaller parts. TOPCONS2 predictions for each part should be concatenated before downstream processing.

Refining TMD boundaries

Refinement of TMD boundaries is performed using von_heijne_scan_functions_v39c_thirdparse.py. For each TOPCONS-predicted transmembrane segment, the script scans windows of length 16–25 residues and identifies the segment with the lowest predicted membrane insertion free energy.

The command used in the manuscript analysis was:

python von_heijne_scan_functions_v39c_thirdparse.py ensembl_homo_sapiens_grch38_15_aug_release_113_translated.fasta topcons_result_ensembl_pep_march_2025_combined.txt topcons -lmin 16 -lmax 25 -o topcons_ensembl_output_march_2025.csv -d dg_topcons_ensembl_output_march_2025.txt

Outputs

The primary output is topcons_ensembl_output_march_2025.csv, which contains refined TMD coordinates and is used as input for the motif-search pipeline.
