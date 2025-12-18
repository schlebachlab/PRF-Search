TMD-slip motif search

This folder contains the core pipeline used to identify TMD-slip motifs and to generate predicted −1 frameshift translation products. A TMD-slip motif is defined as a slippery heptanucleotide sequence with low intrinsic resistance to −1 frameshifting that occurs at a characteristic distance downstream of a transmembrane domain predicted to engage the Sec61 translocon.

The workflow is divided into three stages, each implemented in a subfolder. These stages are typically run in the order listed below, as outputs from earlier steps are required inputs for later steps.

Heptamer Scoring

This stage assigns thermodynamic “friction” scores to candidate heptamer sequences. In the manuscript analysis, all possible heptamers are scored using a Turner nearest-neighbor model. The resulting friction-score table is used to define which sequences are considered slippery in the motif search.

See Heptamer Scoring/README.md.

Transmembrane Domain Analysis

This stage translates Ensembl CDS transcripts into protein sequences, processes TOPCONS2 predictions, and refines transmembrane-domain boundaries using a ΔG-based scan. The output is a standardized CSV file of TMD annotations used by the motif-search pipeline.

See Transmembrane Domain Analysis/README.md.

TMD-slip motif search and frameshifted peptide prediction

This stage performs the transcriptome-wide motif search, associates slippery heptamers with upstream TMDs, applies spacing constraints, and generates both motif tables and predicted frameshift translation products. It also emits digested peptide FASTAs and peptide-to-parent mapping files used in downstream proteomics analyses.

See TMD-slip motif search and frameshifted peptide prediction/README.md.
