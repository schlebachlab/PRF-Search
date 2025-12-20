# PRF-Search
Scripts to search for novel ribosomal frameshifting motifs in human transcripts
TMD-slip motif discovery and frameshift-product generation

This repository contains the analysis code used to identify transmembrane-domain–coupled −1 ribosomal frameshifting (TMD-slip motifs) in the human transcriptome and to generate the predicted frameshift products used in downstream proteomics analyses. The outputs produced by these scripts correspond to those reported in the accompanying manuscript and Supplemental Theory.

The repository is organized around the sequence of analyses performed in the study. Conceptual and theoretical motivation for the scoring models, spacing constraints, and motif definitions are described in the manuscript and supplemental materials. The documentation here focuses on how to reproduce the concrete computational outputs used in the paper.

Conda environment setup
All scripts in this repository were run inside a single conda environment defined by the YAML file located at the repository root.

Create the environment by running:

conda env create -f prf_search_preskon.yml

Activate the environment with:

conda activate prf_search_preskon

Repository layout

The top-level folders correspond to the major stages of the analysis:

TMD-slip motif search
End-to-end pipeline for heptamer friction scoring, transmembrane-domain analysis, motif discovery, and frameshift-product generation.

Statistical Analyses
Scripts used to reproduce the statistical tests and figures reported in the manuscript.

Proteomics Search for Frameshifted Peptides
Utilities for validating spectra assigned to chimeric frameshift products and mapping peptides to parent sequences.

The peptide source utility requires a mapping of peptide to parent produced by the TMD-slip motif search (peptide2parent.tsv.gz). The peptide2parent mapping file is too large for GitHub's single-file Git LFS limits, so it is stored in three parts:

- `Proteomic Search/Peptide Source Tracing Utility/not_trimmed_with_all_zf_v101o_fs_pep_all_msfrag_nme_always_peptide2parent.tsv.gz.part-00`
- `...part-01`
- `...part-02`

To reconstruct the original file:

    bash recombine_peptide2parent.sh

This will create:
Proteomic Search/Peptide Source Tracing Utility/not_trimmed_with_all_zf_v101o_fs_pep_all_msfrag_nme_always_peptide2parent.tsv.gz

Advanced option: regenerate the mapping from source using
TMD-Slip Motif Search/TMD-slip motif search and frameshifted peptide prediction/enhanced_harrington_search_single_ss_per_tmd_v120_fix_motif_table_all.py

Gene Ontology Analysis
Scripts used for GO enrichment and related functional analyses.

To reproduce the primary motif table and the chimeric FASTA files used for mass-spectrometry searches, begin with the TMD-slip motif search folder and follow its README.
