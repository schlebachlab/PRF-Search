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

Gene Ontology Analysis
Scripts used for GO enrichment and related functional analyses.

To reproduce the primary motif table and the chimeric FASTA files used for mass-spectrometry searches, begin with the TMD-slip motif search folder and follow its README.
