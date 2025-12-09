# The script fs_pep_utilities has a number of jobs it can do to query the big table of peptides produced by Harrington motifs and subsequent enzymatic digestion
# The TSV motif_search_multiSS_v99_fs_pep_peptide2parent.tsv contains all predicted fragments produce by digestion, including segments from zero frame (canonical), 
# -1 frame (frameshifted) and peptides spanning the transition from zero frame to -1 frame (transitional)
# These utilities allow us to extract information from this table.
# The table is just a TSV form of the big fasta file motif_search_multiSS_v99_fs_pep_digested_fs_pep.fa
# When I realized these files were unwieldy due to how big they are I decided to write a utility to extract whatever information I thought might be useful
# The TSV and fasta mentioned above are so large because there is a lot of redundant information in them
# This is due to the fact that each peptide can be made in several ways, depending on transcript isoforms, different enzymes, different types of frameshifting, etc
# Therefore the information relevant to each peptide is gathered in a rather long header in the fasta and mapping TSV. The mapping tsv (peptide2parent.tsv) 
# maps each peptide to every one of it's possible headers
# This is the information I want to trim down. We want to ensure that a given peptide can only be made at one of our predicted frameshift sites

# Information in the headers should allow us to trace every peptide back to the TMD-slip motif where it was made, the frameshift type (1-tRNA, 2-tRNA, etc),
# transcript id, gene id, and whether the peptide is from a canonical segment, frameshifted segment, or transitional segment (the novelty= tag in the header has this info)
# We can even use the motif file, motif_search_multiSS_v99_fs_pep_harrington_motifs.csv, to get all the peptides predicted for a motif of interest
# Each motif is a row in that table. Simply copy and paste the header and any row (motif) you are interested in and feed it to this script, using the --motif flag

# Further detail on how to use the utility for these different cases follows.


# Script requires pandas
# conda install pandas
# or
# pip install pandas


# Some peptides originate in the zero reading frame part, some originate in the -1 reading frame segment,
# some have both zero reading frame and -1 reading frame parts (they straddle the point where a ribosomal frameshift
# occurs. If a peptide only has these transitional headers, it is in the set transitional_header; if it only has
# frameshift headers, it is frameshift_only, if it only has zero_frame header tags, it is canonical_only)

# SETS (use --set flag and choose one of these)
#    frameshift_only
#    transitional_only
#    canonical_only
#    transitional_and_frameshift
#    mixture_with_canonical
#    no_zero_frame_header

    # In various older versions of the utility script I called zero_frame canonical; 
    # so if you see --set no_canonical_header or --set canonical_only replace that with --set no_zero_frame_header or --set zero_frame_only

# Some example commands you can run with this utility. 
# --query can be a list of peptides in a .txt file or a fasta file (extention .fa or .fa.gz or .fasta)
# -m should always be the _peptide2parent.tsv.gz
# Tables can be saved as CSV, TSV, or XLSX

# 1) Dump all FS/TS headers for a peptide or set of peptides into a table (TSV/CSV/XLSX); no set filter means that all headers will be printed out
# Headers are expanded in the form of an individual row; therefore many rows can have the same peptide (which is column 1)
# I would always at least filter out the canonical filters with the --set no_canonical_header option
python fs_pep_utilities_v11.py headers -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query drown_fs_peps.txt --out drown_headers_all.tsv

python fs_pep_utilities_v11.py headers -m motif_search_multiSS_v99_fs_pep_peptide2parent_canonical.tsv.gz --query drown_fs_peps.txt --out drown_headers_canonical.tsv

# 1a) Only those headers if the peptide has no canonical (zero reading frame) headers:
 # No query flag means all peptides in the TSV will be written
 # No set flag means all peptides are written no matter their header
python fs_pep_utilities_v11.py headers -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query drown_fs_peps.txt --set no_canonical_header --out drown_no_canonical_headers.tsv


# 2) Count & classify per peptide (with optional summary file)
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query drown_fs_peps.txt --out drown_pep_class_counts_all.tsv --summary-out drown_pep_class_summary_all.txt
Wrote 269 peptides to drown_pep_class_counts_all.tsv

Peptide-set summary
frameshift_only             :     41
transitional_only           :     30
zero_frame_only             :    174
transitional_and_frameshift :      1
mixture_with_zero_frame     :      8
not_found                   :     15



python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_canonical.tsv.gz --query drown_fs_peps.txt --out drown_pep_class_counts_canonical.tsv --summary-out drown_pep_class_summary_canonical.txt

Wrote 269 peptides to drown_pep_class_counts_canonical.tsv

Peptide-set summary
frameshift_only             :     41
transitional_only           :     29
zero_frame_only             :     84
transitional_and_frameshift :      1
mixture_with_zero_frame     :      3
not_found                   :    111
Wrote summary to drown_pep_class_summary_canonical.txt


# 2a) Limit the report to peptides in the “no_canonical_header” category
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query drown_fs_peps.txt --set no_canonical_header --out drown_nc_pep_class_counts.tsv --summary-out drown_nc_pep_class_summary.txt


# 3) Upload a set of rows from the TMD-Slip motif csv (must include header, too) and get a list of all the peptides these motifs produce with accompanying header information; you can also filter by peptide set
python fs_pep_utilities_v11.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --motif motif_data_CCR5_SLC6A4_SLC6A5_KCNQ1_for_testing.csv --out test_motif_peps_CCR5_SLC6A4_SLC6A5_KCNQ1_test3.tsv

# 3a) But only report those motif hits whose peptide is in “transitional_only”:
python fs_pep_utilities_v11.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --motif motif_data_CCR5_SLC6A4_SLC6A5_KCNQ1_for_testing.csv --set transitional_only --out test_motif_ts_only_test3a.tsv


# 4) Build a FASTA of arbitrary peptides, filtered by novelty and/or by set
python fs_pep_utilities_v11.py fasta -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query drown_fs_peps.txt --novelty frameshift transitional --out drown_peps_in_cpk_table.fa
  
# 4a) And restrict that FASTA to only “transitional_and_frameshift” peptides:
python fs_pep_utilities_v11.py fasta -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query drown_fs_peps.txt --novelty frameshift transitional --set no_canonical_header --out drown_fs_or_ts_only_test4a.fa

# 5) Build a non-exploded table summarizing all the information in the exploded table, using a fasta or peptide list
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique.fa.gz --set no_canonical_header --out summary_table_v99_fs_pep_digest_unique_no_zeroframe_test5.tsv

# 6a) Build a FASTA of arbitrary peptides to obtain only canonical fragments
python fs_pep_utilities_v11.py fasta -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique.fa.gz --set canonical_only --dedup --out motif_search_multiSS_v99_fs_pep_digest_unique_canonical_only.fa.gz

# 6b) Build a FASTA of arbitrary peptides to obtain only canonical fragments
python fs_pep_utilities_v11.py fasta -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique.fa.gz --set no_canonical_header --dedup --out motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only.fa.gz

# 7a) Count and classify peptides from non-canonical only fasta
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only.fa.gz --out unique_noncanonical_only_pep_class_summary_test7a_table.tsv --summary-out unique_noncanonical_only_pep_class_summary_test7a.txt
 
 # 7b) Count and classify peptides from canonical only fasta
 python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique_canonical_only.fa.gz --out unique_canonical_only_pep_class_summary_test7b_table.tsv --summary-out unique_canonical_only_pep_class_summary_test7b.txt
 
 # 8) Get the motifs that are unique to alternative isoforms and write them to a table
python fs_pep_utilities_v11.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --set no_zero_frame_header --out motif_fs_peps_v99_all_no_zero_frame_header_only_exploded.tsv

# 9) Write these unique_to_alt peptides to a fasta file with --dedup option
python fs_pep_utilities_v11.py fasta -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query v99_peps_unique_to_alt_isoforms.txt --set no_zero_frame_header --dedup --out motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only_alt_isoform_unique.fa.gz

# 10) Count and classify peptides from motifs that are unique to alternative isoforms (still need to filter these peptides from those in canonical isoforms)
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only_alt_isoform_unique.fa.gz --out non_zero_frame_peptides_from_motifs_unique_to_alt_isoforms_table.tsv --summary-out non_zero_frame_pep_motifs_unique_to_alt_isoforms_summary_out.txt

# 11) Build a summary table of peptides from motifs unique to alt isoforms
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only_alt_isoform_unique.fa.gz --set no_zero_frame_header --out summary_table_v99_fs_ts_peps_from_motifs_unique_to_alt_isoforms.tsv

# 12) Get a fasta of peptides unique to alternative isoforms (with no peptides from other motifs present in the fasta) USES DIFFERENT SCRIPT - fasta_diff.py
python fasta_diff.py motif_search_multiSS_v99_fs_pep_digest_unique_fs_or_ts_only_alt_isoform_unique.fa.gz motif_search_multiSS_v99_fs_pep_digest_unique_all.fa.gz peptides_only_in_motifs_unique_to_alt_isoforms_diff_v99.fa

# 13 Write a table of peptides from motifs unique to alt isoforms
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query peptides_only_in_motifs_unique_to_alt_isoforms_diff_v99.fa --out summary_table_v99_peps_from_motifs_unique_to_alt_isoforms_diff.tsv

# 13a Count and classify peptides from motifs unique to alt isoforms
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query peptides_only_in_motifs_unique_to_alt_isoforms_diff_v99.fa --out peptides_from_motifs_unique_to_alt_isoforms_diff_class_table.tsv --summary-out pep_motifs_unique_to_alt_isoforms_diff_summary_out.txt

# NOTE: tests 8-13 yielded an inaccurate set of peptides. Use 14 below instead. Keeping those commands here in case the other cases are useful

# 14) Build a fasta of peptides that are truly unique to peptides originating motifs unique to alternative isoform transcripts (meaning not in canonical transcripts)
python fs_pep_utilities_v11_dedup.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --filter --fa unique_alt_isoform_peptides_dedup.fa.gz --dedup

# This gave zero motif-driven peptides found. Let's build the fasta of these peptides instead, then manually inspect them with the table
# 14a) Build a fasta as described in number 14 but without eliminating peptides found in other motifs
python fs_pep_utilities_v11_dedup.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --filter --fa unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --dedup
No motif-driven peptides found
# Still no motif-driven peptides found. Let's try version11 without dedup option in case that version of script broke things
# 14b)
python fs_pep_utilities_v11.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --filter --fa unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz

# Still no answer. Okay, let's remove back to v10 and get the list of motifs from that
# 14c) # THIS ONE GAVE AN ANSWER
python fs_pep_utilities_v11_dedup.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --fa unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --dedup

# 14d) do this with --set option
python fs_pep_utilities_v11_dedup.py motif -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --set no_zero_frame_header --motif motifs_v99_L35-55_unique_to_alternative_isoforms.csv --fa unique_to_alt_isoform_peptides_dedup_no_diff_performed_no_zeroframe.fa.gz --dedup

# Okay, that worked. Real unique peptides may actually exist if you reduced missed cleavages parameter.
# Now let's build a summary table of this TSV
# 15) Build a summary table of peptides from motifs that are unique to alternative isoforms (even though they can also be found from other sources apparently)
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --out unique_to_alt_isoform_peptides_dedup_no_diff_performed_summary_table.tsv

# 15-2) Now do this on the version without canonical headers
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed_no_zeroframe.fa.gz --out unique_to_alt_isoform_peptides_dedup_no_diff_no_zeroframe_performed_summary_table.tsv

# 15a) Now run the classify module on it
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --out unique_to_alt_isoform_peptides_dedup_no_diff_performed_class_table.tsv --summary-out pep_motifs_unique_to_alt_isoforms_diff_classify_out.txt

# 15a-2) Now run the classify module on the version without canonical headers
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --out unique_to_alt_isoform_peptides_dedup_no_diff_performed_no_0frame_hdrs_class_table.tsv --summary-out pep_motifs_unique_to_alt_isoforms_diff_no_0frame_hdrs_classify_out.txt

# 16) One last sanity check. Use the script filter_alt_only_peptides.py to get the collection of all peptides with no headers containing 'transcript_canonical_status=Canonical'
python filter_alt_only_peptides.py -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz -o alt_only_peptides.txt
# okay, we are onto something. Here was the output: Wrote 450632 alternative-only peptides to alt_only_peptides.txt

# 17) Now we have a set of peptides. Let's classify them.
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query alt_only_peptides.txt --out alt_only_peps_classified.tsv --summary-out alt_only_peps_classified_summary_output.txt

# Looks like we have some winners!
#Peptide-set summary
#frameshift_only             :  36991
#transitional_only           :  80230
#zero_frame_only             : 330028
#transitional_and_frameshift :    587
#mixture_with_zero_frame     :   2796
#not_found                   :      0
#Wrote summary to alt_only_peps_classified_summary_output.txt

# 18) Make a fasta of this. Nah, these are just all the peptides from alternative transcripts. Not ones from chimeras from motifs that are uniquely present in alternative isoforms.

# 19) Classify your previous fasta generated from motifs that are unique to alt isoforms:
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --out alt_only_motif_peps_classified.tsv --summary-out alt_only_motif_peps_classified_summary_output.txt

Wrote 745,747 peptides to alt_only_motif_peps_classified.tsv

Peptide-set summary
frameshift_only             :  36941
transitional_only           :  88836
zero_frame_only             : 609200
transitional_and_frameshift :    671
mixture_with_zero_frame     :  10099
not_found                   :      0
Wrote summary to alt_only_motif_peps_classified_summary_output.txt

# 20) Repeat classification with just the ones that don't contain zero-frame segments
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed_no_zeroframe.fa.gz --out alt_only_motif_peps_no0framehdrs_classified.tsv --summary-out alt_only_motif_peps_no0framehdrs_classified_summary_output.txt

Wrote 126,448 peptides to alt_only_motif_peps_no0framehdrs_classified.tsv

Peptide-set summary
frameshift_only             :  36941
transitional_only           :  88836
zero_frame_only             :      0
transitional_and_frameshift :    671
mixture_with_zero_frame     :      0
not_found                   :      0
Wrote summary to alt_only_motif_peps_no0framehdrs_classified_summary_output.txt

# 19-2) Classify your previous fasta generated from motifs that are unique to alt isoforms from canonical mapping only:
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_canonical.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed.fa.gz --out alt_only_motif_peps_classified_cn_only.tsv --summary-out alt_only_motif_peps_classified_summary_output_cn_only.txt

20-2) Classify your previous fasta generated from motifs that are unique to alt isoforms so none contain zero frame segments but are from canonical only mapping
python fs_pep_utilities_v11.py classify -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query unique_to_alt_isoform_peptides_dedup_no_diff_performed_no_zeroframe.fa.gz --out alt_only_motif_peps_cn_only_no0framehdrs_classified.tsv --summary-out alt_only_motif_peps_cn_only_no0framehdrs_classified_summary_output.txt

21) Build a summary table of all unique peptides that aren't zero-frame for all peptides
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query motif_search_multiSS_v99_fs_pep_digest_nonzeroframe_only_all.fa.gz --set no_zero_frame_header --out summary_table_v99_fs_pep_digest_all_unique_no_zeroframe.xlsx

22) Build summary table of drown peptides that aren't zero-frame
python fs_pep_utilities_v11.py summarize -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query drown_fs_peps.txt --set no_zero_frame_header --out summary_table_v99_fs_pep_digest_drown_unique_no_zeroframe.xlsx

23) Print out all headers associated with non zero-frame Drown peptides
python fs_pep_utilities_v11.py headers -m motif_search_multiSS_v99_fs_pep_peptide2parent_all.tsv.gz --query drown_fs_peps.txt --set no_zero_frame_header --out all_headers_drown_no_zeroframe_headers.xlsx
