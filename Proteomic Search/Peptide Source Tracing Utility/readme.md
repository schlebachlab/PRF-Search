Peptide Source Tracing Utility

This folder contains fs_pep_utilities_script_updated_enz_rules_v17.py, a command-line utility for tracing peptide-spectrum matches (PSMs) back to the specific chimeric translation products and motif instances they could have come from. In our workflow, MSFragger is run on the full-length chimera FASTA produced by the TMD-slip motif search pipeline. For the peptides that match spectra, this script uses the accompanying peptide2parent.tsv(.gz) mapping (also emitted by the motif search pipeline) to determine whether each peptide is compatible with the 0-frame prefix, a transitional segment spanning the frame junction, or a frameshift-derived segment downstream of the slippery site. The goal is to rapidly identify peptides that are uniquely diagnostic of frameshifting, versus peptides that could also be explained by 0-frame translation.

The core input is a peptide2parent mapping table: one row per peptide, with a serialized header describing its source (transcript/gene identifiers, motif context, novelty class, enzyme/digest coordinates, etc.). The script provides several subcommands. Most commonly, you will (i) classify a peptide list returned by MSFragger into evidence categories and (ii) “explode” the underlying headers for the subset you want to inspect in detail.

Inputs you will need

A peptide-to-parent mapping table from the motif search pipeline (often gzipped), e.g.
not_trimmed_with_all_zf_v101o_fs_pep_all_msfrag_nme_always_peptide2parent.tsv.gz

A “query” file containing peptides returned by MSFragger (one peptide per line works; a TSV with a peptide column also works in our usage). In November 2025 we used, for example:
post_fragpipe_v101o_predicted_frameshift_peptides_21_nov_2025.txt

(Optional) If you want to restrict to particular evidence sets (frameshift-only, transitional-only, etc.), you can filter after classification or use --set on some subcommands. The CLI exposes --set for headers (and related logic internally) and describes the available set names. 

What the utility produces

Depending on the subcommand, outputs include:

A TSV table that assigns each peptide to a high-level evidence category such as frameshift_only, transitional_only, zero_frame_only, mixture_with_zero_frame, etc. (This is the output of classify.) 

A TSV table of exploded header fields for each peptide-source pairing (this is the output of headers). This is what we use to see which transcript/motif instance(s) could generate a peptide. The headers command supports a streaming low-memory mode (--lean) and an optional explicit column list (--cols). 

A FASTA file of peptides (optionally deduplicated) for downstream use (this is the output of fasta). 

A one-row-per-peptide summary table (this is the output of summarize), which is handy when you want counts of the different novelty classes and a compact list of associated transcripts/genes/motif metadata. 

The primary script in this folder is fs_pep_utilities_script_updated_enz_rules_v17.py.

Core use case: classifying MSFragger peptide hits

In the manuscript workflow, this utility is used after MSFragger has searched spectra against a FASTA of full-length chimeric translation products. For peptides that match spectra, the collaborator provides a list of peptide sequences. These peptides are then classified to determine whether they provide evidence for frameshifting.

The classifier is run as follows:

python fs_pep_utilities_script_updated_enz_rules_v17.py classify -m not_trimmed_with_all_zf_v101o_fs_pep_all_msfrag_nme_always_peptide2parent.tsv.gz --query post_fragpipe_v101o_predicted_frameshift_peptides_21_nov_2025.txt --out post_fragpipe_v101o_predicted_frameshift_peptides_21_nov_2025_peps_classified.tsv --summary-out post_fragpipe_v101o_predicted_frameshift_peptides_21_nov_2025_classified_summary.tsv --not-found-table post_fragpipe_v101o_predicted_frameshift_peptides_21_nov_2025_classified_not_found.tsv

Here, the mapping file (-m) encodes all possible peptide sources from the chimeric translation products, and the query file (--query) contains the peptide sequences identified by MSFragger.

Outputs

The primary output is a TSV file (--out) containing one row per queried peptide, annotated with counts of compatible parent sequences and an overall classification label. Classification categories include:

frameshift_only

transitional_only

zero_frame_only

transitional_and_frameshift

zero_frame_with_frameshift_stop_only

mixture_with_zero_frame

not_found

A summary file (--summary-out) is also written, reporting the total number of peptides assigned to each category. For the 21 November 2025 analysis, the summary was:

frameshift_only : 14
transitional_only : 17
zero_frame_only : 462
transitional_and_frameshift : 0
zero_frame_with_frameshift_stop_only: 7
mixture_with_zero_frame : 13
not_found : 0

An additional table (--not-found-table) lists any peptides that could not be matched to entries in the mapping file, although this is typically empty when the mapping and FASTA inputs are consistent.

Interpretation

Peptides classified as frameshift_only or transitional_only cannot be explained by conventional 0-frame translation and therefore constitute direct proteomic evidence for frameshifting at TMD-slip motifs. These classifications are used in the manuscript to support the existence of frameshift-derived protein products.
