#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to search transcriptome databases (Ensembl CDS) for Harrington motifs
This fork creates fasta files of 3 types frameshift products at slippery sequences
for each identified motif as well as peptide fragments of these after treatment
with several proteases
"""
import RNA
import re
import textwrap
import csv
import numpy as np
import pandas as pd
import random
import sys
import pickle
import gc
from collections import Counter, defaultdict
from pyteomics import parser
import multiprocessing as mp
from scipy import stats
from scipy.stats import mannwhitneyu, binomtest, ks_2samp, kstwobign, zscore, chisquare
from scipy.stats import norm
from itertools import groupby
import matplotlib.pyplot as plt
import frameshift_routines_v10 as frameshift
import time
import configparser
from Bio import Align
from Bio.Align import substitution_matrices
import difflib
import gzip
from typing import Tuple

from dataclasses import dataclass
@dataclass
class DigestOptions:
    # enzyme behavior
    gluc_mode: str = "E_or_D"           # {"E_only","E_or_D"}
    chymo_specificity: str = "core"     # {"core","low"}  core=FYW, low=FYW+L
    trypsin_mode: str = "trypsin"       # {"trypsin","stricttrypsin"}  (both respect K/R|P exception)
    # search behavior
    semi_enzymatic: bool = True         # allow one non-enzymatic terminus
    missed: int = 2
    min_len: int = 6
    max_len: int = 50
    # proteoform behavior
    nme_mode: str = "bio"               # {"off","bio","always"}

# Global variables
friction_data = None
scoring_sys = None



# Functions

def parse_blastp_output(blastp_results_path):
    """
    Parses the blastp output (TSV) and returns a dictionary where keys are Ensembl transcript ID and values are
    a list of matches, each with Uniprot Accession ID and similarity score.
    """
    blastp_matches = {}
    with open(blastp_results_path, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            ensembl_transcript_id = row[0]
            uniprot_accession_id = row[1].split('|')[1]  # Extract the Uniprot ID from the full string
            similarity_score = float(row[2])

            # Add the match to the dictionary
            if ensembl_transcript_id not in blastp_matches:
                blastp_matches[ensembl_transcript_id] = []
            blastp_matches[ensembl_transcript_id].append((uniprot_accession_id, similarity_score))
    
    return blastp_matches


def import_tsv(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.replace('\n', '').split('\t'))
    file_obj.close()
    return data_raw


def parse_ensembl_uniprot_canon(path_ensembl_uniprot_tsv):
    '''
    'Gene stable ID',               0
    'Gene stable ID version',       1
    'Transcript stable ID',         2
    'Transcript stable ID version', 3
    'UniProtKB/Swiss-Prot ID',      4
    'Gene name',                    5
    'Transcript type',              6
    'Ensembl Canonical'             7
    '''

    ensembl_features_dict = {}

    ensembl_uniprot_features_raw = import_tsv(path_ensembl_uniprot_tsv)
    ensembl_uniprot_features_list = ensembl_uniprot_features_raw[1:]

    for entry in ensembl_uniprot_features_list:
        stable_transcript_id = entry[2]
        stable_gene_id = entry[0]
        transcript_type = entry[6]
        transcript_id_version = entry[3]

        if entry[4] != '':
            uniprot_id = entry[4]
        else:
            uniprot_id = 'No_uniprot_annotation'

        if entry[5] != '':
            gene_name = entry[5]
        else:
            gene_name = 'No_gene_name_annotation'

        if entry[7] == '1':
            canonical_status = 'Canonical'
        else:
            canonical_status = 'Alternative'

        dict_entry = [transcript_id_version, stable_gene_id, transcript_type, uniprot_id, gene_name, canonical_status]
        # transcript_id_version = ensembl_features_dict[stable_transcript_id][0][0]
        # stable_gene_id = ensembl_features_dict[stable_transcript_id][0][1]
        # transcript_type = ensembl_features_dict[stable_transcript_id][0][2]
        # uniprot_id = ensembl_features_dict[stable_transcript_id][0][3]
        # gene_name = ensembl_features_dict[stable_transcript_id][0][4]
        # canonical_status = ensembl_features_dict[stable_transcript_id][0][5]
        # [transcript_id_version, stable_gene_id, transcript_type, uniprot_id, gene_name, canonical_status] = ensembl_features_dict[stable_transcript_id][0]

        if stable_transcript_id not in ensembl_features_dict:
            ensembl_features_dict[stable_transcript_id] = []
        ensembl_features_dict[stable_transcript_id].append(dict_entry)

    return ensembl_features_dict



def fasta_transcript_filter(fasta_list, ensembl_features_dict):
    # This is the indexing of the output
    # output_line = [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len(sequence_str), sequence_str, canonical_status, uniprot_id, transcript_type, alias_transcript_ids]
    #                    0                 1                2           3                 4                 5             6               7               8               9
    print('Running Ensembl transcript filter...')

    unique_genes = []
    initial_labeled_list = []
    labeled_list = []


    # Identify unique genes
    for entry in fasta_list:
        current_gene = entry[0]
        if current_gene not in unique_genes:
            unique_genes.append(current_gene)
    
    # Apply features from ensembl feature dict to each entry
    for entry in fasta_list:
        gene_id_version = entry[0]
        gene_id_stable = gene_id_version.split('.')[0]
        fasta_gene_name = entry[2]
        transcript_id_version = entry[1]
        transcript_id_stable = transcript_id_version.split('.')[0]
        if transcript_id_stable not in ensembl_features_dict:
            ensembl_features_dict[transcript_id_stable] = [[transcript_id_version, gene_id_stable, 'transcript_type_unknown' ,'No_uniprot_annotation' , fasta_gene_name, 'canonical_status_unknown']]
        entry_ensembl_features = ensembl_features_dict[transcript_id_stable]
        if len(entry_ensembl_features) == 1:
            [transcript_id_full, stable_gene_id, transcript_type, uniprot_id, gene_name, canonical_status] = entry_ensembl_features[0]
            entry.append(canonical_status)
            entry.append(uniprot_id)
            entry.append(transcript_type)
        else:
            [transcript_id_full, stable_gene_id, transcript_type, uniprot_id_first, gene_name, canonical_status] = entry_ensembl_features[0]
            entry.append(canonical_status)
            sub_entry_uniprot_ids = []
            for sub_entry in entry_ensembl_features:
                sub_entry_uniprot_id = sub_entry[3]
                if sub_entry_uniprot_id not in sub_entry_uniprot_ids:
                    sub_entry_uniprot_ids.append(sub_entry[3])
            if len(sub_entry_uniprot_ids) > 1:
                uniprot_id = '; '.join(sub_entry_uniprot_ids)
            else:
                uniprot_id = uniprot_id_first
            entry.append(uniprot_id)
            entry.append(transcript_type)

        initial_labeled_list.append(entry)

    # Create bins for transcripts belonging to the same gene
    binned_list = [[] for _ in range(len(unique_genes))]
    for idx, gene_id in enumerate(unique_genes):
        for entry in initial_labeled_list:
            if entry[0] == gene_id:
                binned_list[idx].append(entry)

    # Check gene bins to ensure that there is 1 canonical transcript per ensembl gene
    genes_with_no_canonical_transcripts = []
    genes_with_one_canonical_transcript = []
    genes_with_multiple_canonical_transcripts = []

    for gene_bin in binned_list:
        num_canonical = 0
        for transcript_entry in gene_bin:
            canonical_status = transcript_entry[6]
            if canonical_status == 'Canonical':
                num_canonical += 1
        if num_canonical == 1:
            genes_with_one_canonical_transcript.append(gene_bin)
        if num_canonical == 0:
            genes_with_no_canonical_transcripts.append(gene_bin)
        if num_canonical > 1:
            genes_with_multiple_canonical_transcripts.append(gene_bin)




    for gene_bin in genes_with_one_canonical_transcript:
        for entry in gene_bin:
            labeled_list.append(entry)

 
    dedup_list = []

    # Routines for de-duplication
    # Create bins for transcripts belonging to the same gene.
    binned_dict = {}  # key: gene ID, value: list of transcript entries
    for entry in labeled_list:
        gene_id = entry[0]
        if gene_id not in binned_dict:
            binned_dict[gene_id] = []
        binned_dict[gene_id].append(entry)


    # Now, for each gene, remove redundancy based on CDS sequence (entry[5])
    for gene_id, gene_bin in binned_dict.items():
        # Create a dictionary mapping CDS sequence to the list of transcript entries that have that sequence.
        seq_to_entries = {}
        for entry in gene_bin:
            seq = entry[5]
            if seq not in seq_to_entries:
                seq_to_entries[seq] = []
            seq_to_entries[seq].append(entry)


        # For each unique sequence, choose a representative and record aliases
        for seq, entries in seq_to_entries.items():
            if len(entries) == 1: # Here entries is a list of entries that is the value to unique sequences, which are keys
                representative_entry = entries[0]
                # Add an alias field; if no duiplicates assign 'no_alias'
                representative_entry.append("no_alias")
                dedup_list.append(representative_entry)
            else:
                # There are duplicate transcripts with identical CDS
                # If any transcript is canonical, choose that transcript as representative
                canonicals = [entry for entry in entries if entry[6] == 'Canonical']
                if canonicals:
                    representative_entry = canonicals[0] # Choose the first canonical encountered (should only be one in labeled_list)
                else:
                    representative_entry = entries[0] # Choose the first one if none are canonical
                # Build the alias string from the transcript IDs (stable IDs with version) of the others
                representative_transcript_id = representative_entry[1]
                alias_transcript_ids = []
                for entry in entries:
                    current_transcript_id = entry[1]
                    if current_transcript_id != representative_transcript_id:
                        alias_transcript_ids.append(current_transcript_id)
                if alias_transcript_ids:
                    alias_str = '; '.join(alias_transcript_ids)
                else:
                    alias_str = 'no_alias'
                representative_entry.append(alias_str)
                dedup_list.append(representative_entry)


    output_list = []
    # Only use this loop if you are using the set genes_with_one_canonical_transcript
    # This goes through the output list and changes any alternative isoforms that were labeled unknown to 'Alternative'
    # We can do this if we know that each gene has one canonical transcript
    #for entry in labeled_list: # Un-comment this line and recomment the one below if you want to use labeled list
    for entry in dedup_list: # Comment out this line and un-comment the one above if you want to use labeled list
        canonical_status = entry[6]
        if canonical_status == 'canonical_status_unknown':
            entry[6] = 'Alternative'
        output_list.append(entry)

    num_canonical_transcripts = 0
    num_alternative_transcripts = 0
    num_canonical_unknown_transcripts = 0

    #for entry in initial_labeled_list: # If you want diagnostics on how many outliers there are uncomment this and comment the line below
    #for entry in labeled_list:
    for entry in output_list:
        canonical_status = entry[6]
        if canonical_status == 'Canonical':
            num_canonical_transcripts += 1
        if canonical_status == 'Alternative':
            num_alternative_transcripts += 1
        if canonical_status == 'canonical_status_unknown':
            num_canonical_unknown_transcripts += 1

    # Comment out this line if you don't want to learn how many duplicates were eliminated
    num_duplicates_eliminated = len(labeled_list) - len(dedup_list)
    print('Number of duplicate sequences filtered: ', num_duplicates_eliminated)
    print('Number of transcripts obtained from Ensembl after filters and de-duplication: ', len(output_list))
    print('Number of canonical transcripts identified: ', num_canonical_transcripts)
    print('Number of alternative transcripts identified: ', num_alternative_transcripts)

    # diagnostic variables returned in this commented out return line
    #return initial_labeled_list, binned_list, genes_with_no_canonical_transcripts, genes_with_one_canonical_transcript, genes_with_multiple_canonical_transcripts, num_canonical_transcripts, num_alternative_transcripts, num_canonical_unknown_transcripts

    return output_list, num_canonical_transcripts, num_alternative_transcripts



def bin_transcripts_by_gene_id(filtered_list):
    """
    Bins Ensembl transcripts by their shared Ensembl Gene ID.
    
    Parameters:
    filtered_list (list): List of parsed and filtered Ensembl transcripts.
    
    Returns:
    dict: Dictionary where keys are Ensembl Gene IDs and values are lists of transcripts associated with that Gene ID.
    """
    binned_dict = {}
    
    # Iterate through the filtered_list and bin by ensembl_gene_id (which is the first item in each entry)
    for entry in filtered_list:
        ensembl_gene_id = entry[0]  # Assuming the Ensembl Gene ID is at index 0
        
        if ensembl_gene_id not in binned_dict:
            binned_dict[ensembl_gene_id] = []
        
        # Append the transcript entry to the list of transcripts for this gene
        binned_dict[ensembl_gene_id].append(entry)
    
    return binned_dict


def match_ensembl_to_uniprot(blastp_results, ensembl_gene_id_dict):
    """
    Matches Ensembl transcripts to Uniprot sequences using the parsed blastp output.
    
    Arguments:
    - filtered_list: List of Ensembl transcripts, labeled as canonical or alternative.
    - blastp_results: Dictionary of blastp results, mapping Ensembl transcript IDs to Uniprot matches.
    - ensembl_gene_id_dict: Dictionary mapping Ensembl Gene IDs to their corresponding transcripts.
    
    Returns:
    - ensembl_to_uniprot: A dictionary mapping Ensembl transcript IDs to Uniprot Accession IDs and similarity scores.
    """
    ensembl_to_uniprot = {}

    for ensembl_gene_id, transcripts in ensembl_gene_id_dict.items():
        # Keep track of the best Uniprot match (highest similarity score)
        best_match = None
        best_similarity_score = -1
        best_uniprot_id = None
        
        # First pass: identify the best Uniprot match for the canonical transcript
        for transcript in transcripts:
            transcript_id = transcript[1]  # Ensembl transcript ID

            if transcript_id in blastp_results:
                for uniprot_id, similarity_score in blastp_results[transcript_id]:
                    if similarity_score > best_similarity_score:
                        best_similarity_score = similarity_score
                        best_uniprot_id = uniprot_id
                        best_match = transcript  # Store the transcript with the best match
        
        # If a match is found, update all transcripts in this gene with the Uniprot ID of the canonical transcript
        if best_uniprot_id:
            for transcript in transcripts:
                transcript_id = transcript[1]  # Ensembl transcript ID
                # Check if transcript has an individual match
                transcript_matches = blastp_results.get(transcript_id)
                if transcript_matches:
                    for uniprot_id, similarity_score in transcript_matches:
                        if uniprot_id == best_uniprot_id:
                            ensembl_to_uniprot[transcript_id] = [uniprot_id, similarity_score]
                            break
                    else:
                        ensembl_to_uniprot[transcript_id] = [best_uniprot_id, "No match to this alternative transcript"]
                else:
                    # No match found for this transcript
                    ensembl_to_uniprot[transcript_id] = [best_uniprot_id, "No match to this alternative transcript"]
        else:
            # No match found for any transcript in this gene
            for transcript in transcripts:
                transcript_id = transcript[1]
                ensembl_to_uniprot[transcript_id] = ["No Uniprot Accession ID found", "No similarity score with any Uniprot sequence found"]

    return ensembl_to_uniprot



def process_transcripts_with_uniprot(filtered_list, blastp_results_path, tmd_data, outname):
    """
    Processes Ensembl transcripts to associate them with Uniprot sequences based on blastp results.

    Parameters:
    - filtered_list (list): The list of parsed and filtered Ensembl transcripts.
    - blastp_results_path (str): Path to the TSV file with blastp results.

    Returns:
    dict: A dictionary mapping Ensembl transcript IDs to Uniprot metadata, including Uniprot accession and similarity score.
    """
    # Step 1: Parse the blastp results
    blastp_results = parse_blastp_output(blastp_results_path)

    # Step 2a: Filter transcripts that lack a TMD
    # Skip this for now - will filter by KW-0472 later after harrington motif list has been made
    #tmd_filtered_list = filter_transcripts_with_tmd(filtered_list, tmd_data)

    # Step 2b: Bin transcripts by Ensembl Gene ID
    ensembl_gene_id_dict = bin_transcripts_by_gene_id(filtered_list)

    # Step 3: Match Ensembl transcripts to Uniprot sequences
    ensembl_uniprot_dict = match_ensembl_to_uniprot(blastp_results, ensembl_gene_id_dict)

    output_file = outname+'_ensembl_uniprot_dict_dataframe.xlsx'

    # Step 4: Write the results to an Excel file
    output_ensembl_to_uniprot_excel(ensembl_uniprot_dict, filtered_list, output_file)

    return ensembl_uniprot_dict


def output_ensembl_to_uniprot_excel(ensembl_uniprot_dict, filtered_list, output_file):
    """
    Outputs the mapping between Ensembl transcripts and Uniprot sequences to an Excel file.

    Parameters:
    - ensembl_uniprot_dict (dict): A dictionary mapping Ensembl transcript IDs to Uniprot accession IDs and similarity scores.
    - filtered_list (list): The list of parsed and filtered Ensembl transcripts.
    - output_file (str): Path to save the output Excel file.
    """
    data = []
    
    for entry in filtered_list:
        ensembl_gene_id = entry[0]  # Assuming Ensembl Gene ID is at index 0
        ensembl_transcript_id = entry[1]  # Assuming Ensembl Transcript ID is at index 1
        canonical_status = entry[6]  # Assuming canonical status is at index 6

        # Retrieve the Uniprot data
        uniprot_data = ensembl_uniprot_dict.get(ensembl_transcript_id, ["No Uniprot Accession ID found", "No similarity score found"])
        uniprot_accession_id = uniprot_data[0]
        similarity_score = uniprot_data[1]

        # Append the row to the data list
        data.append([ensembl_transcript_id, ensembl_gene_id, canonical_status, uniprot_accession_id, similarity_score])

    # Convert to a DataFrame
    df = pd.DataFrame(data, columns=['Ensembl Transcript ID', 'Ensembl Gene ID', 'Canonical Status', 'Uniprot Accession ID', 'Similarity Score'])

    # Output the DataFrame to Excel
    df.to_excel(output_file, index=False)
    print(f"Ensembl to Uniprot mapping saved to {output_file}")


def append_uniprot_accession_to_motif(harrington_motifs, ensembl_uniprot_dict):
    """
    Appends Uniprot accession and similarity score to the Harrington motif data.

    Parameters:
    harrington_motifs (list): List of Harrington motif entries, where each entry has the Ensembl transcript ID as the first element.
    ensembl_uniprot_dict (dict): Dictionary mapping Ensembl transcript IDs to a list containing Uniprot accession and similarity score.

    Returns:
    list: Updated Harrington motifs list with Uniprot accessions appended to each entry.
    """
    updated_harrington_motifs = []

    for motif_entry in harrington_motifs:
        ensembl_transcript_id = motif_entry[0]  # Assuming the first column is Ensembl transcript ID
        uniprot_info = ensembl_uniprot_dict.get(ensembl_transcript_id)

        if uniprot_info:
            uniprot_id = uniprot_info[0] if len(uniprot_info) > 0 else 'No Uniprot sequence identified'
            similarity_score = uniprot_info[1] if len(uniprot_info) > 1 else 'N/A'
        else:
            uniprot_id = 'No Uniprot sequence identified'
            similarity_score = 'N/A'

        # Append Uniprot accession and similarity score to the motif entry
        updated_entry = motif_entry + [uniprot_id, similarity_score]
        updated_harrington_motifs.append(updated_entry)

    return updated_harrington_motifs


def genetic_code_func():
    genetic_code = """
        TTT F      CTT L      ATT I      GTT V
        TTC F      CTC L      ATC I      GTC V
        TTA L      CTA L      ATA I      GTA V
        TTG L      CTG L      ATG M      GTG V
        TCT S      CCT P      ACT T      GCT A
        TCC S      CCC P      ACC T      GCC A
        TCA S      CCA P      ACA T      GCA A
        TCG S      CCG P      ACG T      GCG A
        TAT Y      CAT H      AAT N      GAT D
        TAC Y      CAC H      AAC N      GAC D
        TAA *      CAA Q      AAA K      GAA E
        TAG *      CAG Q      AAG K      GAG E
        TGT C      CGT R      AGT S      GGT G
        TGC C      CGC R      AGC S      GGC G
        TGA *      CGA R      AGA R      GGA G
        TGG W      CGG R      AGG R      GGG G
    """
    codon_finder = re.compile(r'[ATCG]{3}')
    amino_acid_finder = re.compile(r'\ \w{1}[\ |\n]|\*')
    codon_list = codon_finder.findall(genetic_code)
    amino_acid_list = [x.strip() for x in amino_acid_finder.findall(genetic_code)]
    genetic_code_dict = {}
    i = 0
    while i < len(codon_list):
        genetic_code_dict[codon_list[i]] = amino_acid_list[i]
        i += 1
    return genetic_code_dict


def ribosome(domain_codons, genetic_code_dict):
    domain_aminoacids = []
    for codon in domain_codons:
        amino_acid = genetic_code_dict[codon]
        domain_aminoacids.append(amino_acid)
    return domain_aminoacids


def get_codon_table():
    codon_table = {
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'C': ('TGT', 'TGC'),
        'D': ('GAT', 'GAC'),
        'E': ('GAA', 'GAG'),
        'F': ('TTT', 'TTC'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'I': ('ATT', 'ATC', 'ATA'),
        'H': ('CAT', 'CAC'),
        'K': ('AAA', 'AAG'),
        'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'M': ('ATG',),
        'N': ('AAT', 'AAC'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'Q': ('CAA', 'CAG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG'),
        'W': ('TGG',),
        'Y': ('TAT', 'TAC'),
        '*': ('TAA', 'TAG', 'TGA'),
    }
    return codon_table


def get_codon_freq_table():
    codon_freq_table = {
    'A': (['GCT', 0.26], ['GCC', 0.40], ['GCA', 0.23], ['GCG', 0.11]),
    'C': (['TGT', 0.45], ['TGC', 0.55]),
    'D': (['GAT', 0.46], ['GAC', 0.54]),
    'E': (['GAA', 0.42], ['GAG', 0.58]),
    'F': (['TTT', 0.45], ['TTC', 0.55]),
    'G': (['GGT', 0.16], ['GGC', 0.34], ['GGA', 0.25], ['GGG', 0.25]),
    'I': (['ATT', 0.36], ['ATC', 0.48], ['ATA', 0.16]),
    'H': (['CAT', 0.41], ['CAC', 0.59]),
    'K': (['AAA', 0.42], ['AAG', 0.58]),
    'L': (['TTA', 0.07], ['TTG', 0.13], ['CTT', 0.13], ['CTC', 0.20], ['CTA', 0.07], ['CTG', 0.41]),
    'M': (['ATG', 1.0],),
    'N': (['AAT', 0.46], ['AAC', 0.54]),
    'P': (['CCT', 0.28], ['CCC', 0.33], ['CCA', 0.27], ['CCG', 0.11]),
    'Q': (['CAA', 0.25], ['CAG', 0.75]),
    'R': (['CGT', 0.08], ['CGC', 0.19], ['CGA', 0.11], ['CGG', 0.21], ['AGA', 0.20], ['AGG', 0.20]),
    'S': (['TCT', 0.18], ['TCC', 0.22], ['TCA', 0.15], ['TCG', 0.06], ['AGT', 0.15], ['AGC', 0.24]),
    'T': (['ACT', 0.24], ['ACC', 0.36], ['ACA', 0.28], ['ACG', 0.12]),
    'V': (['GTT', 0.18], ['GTC', 0.24], ['GTA', 0.11], ['GTG', 0.47]),
    'W': (['TGG', 1.0],),
    'Y': (['TAT', 0.43], ['TAC', 0.57]),
    '*': (['TAA', 0.28], ['TAG', 0.20], ['TGA', 0.52]),
    }
    return codon_freq_table


def get_GC(sequence):
    num_GC = 0
    seqlen = len(sequence)
    for base in sequence:
        if base == 'G' or base == 'C':
            num_GC += 1
    percent_GC = (num_GC / seqlen) * 100
    return percent_GC


def calc_avg_freq(domain_codons):
    codon_freq_table = get_codon_freq_table()
    codon_freqs_dict = {}
    for residue, codons in codon_freq_table.items():
        for codon, freq in codons:
            codon_freqs_dict[codon] = freq
    domain_freqs = []
    for position in domain_codons:
        domain_freqs.append(codon_freqs_dict[position])
    len_domain = len(domain_codons)
    domain_freqs_sum = sum(domain_freqs)
    domain_freqs_avg = domain_freqs_sum / len_domain
    return domain_freqs_avg


def load_fasta_data(path_fasta):
    with open(path_fasta, "r") as file_handle:
        faiter = (x[1] for x in groupby(file_handle, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = header.__next__()[0:].strip()
            seq = ''.join(s.strip() for s in faiter.__next__())
            yield (header_str, seq)




def parse_fasta(path_fasta):
    fasta_generator = load_fasta_data(path_fasta)
    fasta_list = []
    for header, sequence in fasta_generator:
        header_list = header.split()
        ensembl_gene = header_list[3][5:]
        ensembl_transcript = header_list[0][1:]
        try:
            gene_name = header_list[6][12:]
        except IndexError:
            gene_name = ensembl_gene
        locus_coord_str = header_list[2][11:]
        sequence_str = ''.join(sequence)
        output_line = [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len(sequence_str), sequence_str]
        len_seq = len(sequence_str)
        # Basic checks on quality of transcript:
        pass_conditions = [(len_seq % 3 == 0),
                           (sequence_str[0:3] == 'ATG'),
                           (sequence_str[-3:] in ('TAG', 'TAA', 'TGA')),
                           ('N' not in sequence_str)]
        if all(pass_conditions):
            fasta_list.append(output_line)
    return fasta_list



def stop_codon_filter(cds_fasta_list):
    print('Filtering transcripts with internal stop codons...')
    stop_filtered_list = []
    genetic_code_dict = genetic_code_func()
    for transcript_data in cds_fasta_list:
        sequence_raw = transcript_data[5]
        transcript_seq = re.sub("[^a-zA-Z]", "", sequence_raw)
        seq_codons = textwrap.wrap(transcript_seq, 3)
        seq_protein = ribosome(seq_codons, genetic_code_dict)
        seq_protein_str = ''.join(seq_protein[0:-1])
        if '*' not in seq_protein_str:
            stop_filtered_list.append(transcript_data)
        else:
            continue
    return stop_filtered_list



def fasta_transcript_filter_canonical_only(fasta_list):
    unique_genes = []
    filtered_list = []
    for entry in fasta_list:
        current_gene = entry[2]
        if current_gene not in unique_genes:
            unique_genes.append(current_gene)
    binned_list = []
    for idx in range(len(unique_genes)):
        binned_list.append([])
    for idx in range(len(unique_genes)):
        gene_name = unique_genes[idx]
        for entry in fasta_list:
            if entry[2] == gene_name:
                binned_list[idx].append(entry)
            else:
                continue
    for gene_bin in binned_list:
        len_longest_transcript = 0
        idx_longest_transcript = 0
        for idx in range(len(gene_bin)):
            len_transcript = gene_bin[idx][4]
            if len_transcript > len_longest_transcript:
                len_longest_transcript = len_transcript
                idx_longest_transcript = idx
        filtered_list.append(gene_bin[idx_longest_transcript])
    return filtered_list



def import_text_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    seqs_str_list = []
    for line in data_raw:
        item = line[0]
        seqs_str_list.append(item)
    file_obj.close()
    return seqs_str_list

def import_raw(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    return data_raw


# Existing function to read friction data
def parse_friction_csv(path_friction_csv):
    friction_data = {}
    with open(path_friction_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            heptamer = row['sequence']
            friction_score = float(row['friction_score'])
            friction_data[heptamer] = friction_score
    return friction_data

# New function to predict secondary structure
def predict_secondary_structure(sequence):
    rna_sequence = sequence.replace('T', 'U')
    structure, mfe = RNA.fold(rna_sequence)
    return structure, mfe


def detect_pseudoknot(structure):
    return any(char in structure for char in ['[', ']', '{', '}', '<', '>'])


def detect_stem_loop(structure, min_stem_length=3, max_loop_size=8):
    stack = []
    stems = []

    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                stems.append((start, i))

    # Filter stems by length and check for loops
    for stem in stems:
        stem_length = stem[1] - stem[0] + 1
        if stem_length >= 2 * min_stem_length:  # Each stem has at least min_stem_length pairs
            loop_size = min(structure[stem[0]+1:stem[1]].count('.'), max_loop_size)
            if loop_size > 0:
                return True  # Found a valid stem-loop
    return False


def filter_transcripts_with_tmd(stop_filtered_list, tmd_data):
    """
    Filters stop_filtered_list to only include transcripts that have TMD data in tmd_data.

    Parameters:
    stop_filtered_list (list): List of Ensembl transcripts that passed the stop codon filter.
    tmd_data (dict): Dictionary containing TMD information, with Ensembl transcript IDs as keys.

    Returns:
    list: Filtered list of Ensembl transcripts that have corresponding TMD data.
    """
    filtered_list_with_tmd = []

    for transcript in stop_filtered_list:
        ensembl_transcript_id = transcript[1]  # Assuming the second item in each entry is the Ensembl transcript ID
        
        # Only include the transcript if it has corresponding TMD data in tmd_data
        if ensembl_transcript_id in tmd_data:
            filtered_list_with_tmd.append(transcript)
    
    return filtered_list_with_tmd


def get_canonical_count(filtered_list, tmd_data, harrington_motifs, outname):

    # Total number of transcripts
    total_canonical_count = 0
    total_alt_spliced_count = 0
    total_count = len(filtered_list)

    # Total number of transripts with TMDs
    canonical_count_with_tmds = 0
    alt_spliced_count_with_tmds = 0
    total_count_with_tmds = 0

    # Number of TMDs in topcons data
    unique_tmd_entries = []
    transcripts_with_tmds = []
    for transcript_id, topcons_info in tmd_data.items():
        if transcript_id not in transcripts_with_tmds:
            transcripts_with_tmds.append(transcript_id)

        for idx, tmd_info in enumerate(topcons_info):
            tmd_seq = tmd_info['adjusted_tmd_seq']  # Extract tmd_type
            tm_start = tmd_info['adjusted_tmd_start_position']
            tm_end = tmd_info['adjusted_tmd_end_position']
            tmd_entry = [transcript_id, tmd_seq, tm_start, tm_end]
            unique_tmd_entries.append(tmd_entry)


    num_unique_tmds = len(unique_tmd_entries)
    num_transcripts_with_tmds = len(transcripts_with_tmds)

    # Total number of unique genes in CDS database
    unique_gene_ids = []

    # Total number of unique genes in CDS database with TMDs
    unique_gene_ids_with_tmds = []

    # Total number of unique HGNC symbols in CDS database
    unique_hgnc_symbols = []

    # Total number of unique HGNC symbols in CDS database with TMDs
    unique_hgnc_symbols_with_tmds = []

    unique_transcripts_with_tmds = []

    for entry in filtered_list:
        gene_id = entry[0]
        hgnc_symbol = entry[2]
        transcript_id = entry[1]
        canon_state = entry[6]
        if gene_id not in unique_gene_ids:
            unique_gene_ids.append(gene_id)

        if hgnc_symbol not in unique_hgnc_symbols:
            unique_hgnc_symbols.append(hgnc_symbol)

        if canon_state == 'Canonical':
            total_canonical_count += 1
        else:
            total_alt_spliced_count += 1

        if transcript_id in tmd_data:
            unique_transcripts_with_tmds.append(transcript_id)
            total_count_with_tmds += 1

            if gene_id not in unique_gene_ids_with_tmds:
                unique_gene_ids_with_tmds.append(gene_id)
            
            if hgnc_symbol not in unique_hgnc_symbols_with_tmds:
                unique_hgnc_symbols_with_tmds.append(hgnc_symbol)

            if canon_state == 'Canonical':
                canonical_count_with_tmds += 1
            else:
                alt_spliced_count_with_tmds += 1
    
    tmd_entries_in_filtered_data = []
    for tmd_entry in unique_tmd_entries:
        transcript_id = tmd_entry[0]
        if transcript_id in unique_transcripts_with_tmds:
            tmd_entries_in_filtered_data.append(tmd_entry)

    num_tmd_entries_in_filtered_data = len(tmd_entries_in_filtered_data)

    # Total number of unique genes
    num_unique_genes = len(unique_gene_ids)

    # Total number of unique genes with TMDs
    num_unique_genes_with_tmds = len(unique_gene_ids_with_tmds)

    # Total number of unique HGNC symbols
    num_unique_hgnc_symbols = len(unique_hgnc_symbols)

    # Total number of unique HGNC symbols with TMDs
    num_unique_hgnc_symbols_with_tmds = len(unique_hgnc_symbols_with_tmds)

    # Number of Harrington motifs identified
    num_harrington_motifs = len(harrington_motifs)

    # Number of Harrington motifs in alternative isoforms
    alt_motifs = 0

    # Number of Harrington motifs in canonical isoforms
    canonical_motifs = 0

    for entry in harrington_motifs:
        canon_status = entry[26]
        if canon_status == 'Alternative':
            alt_motifs += 1
        if canon_status == 'Canonical':
            canonical_motifs += 1

    logfile = open(outname+'_logfile.txt', 'w')
    
    print(f"Total number of transcripts in Ensembl CDS database: {total_count}", file=logfile)
    print(f"Number of canonical transcripts: {total_canonical_count}", file=logfile)
    print(f"Number of alternate isoform transcripts: {total_alt_spliced_count}", file=logfile)
    print(f"Total number of transcripts in Ensembl with TMDs (filtered): {total_count_with_tmds}", file=logfile)
    print(f"Number of canonical transcripts with TMDs: {canonical_count_with_tmds}", file=logfile)
    print(f"Number of alternate isoform transcripts with TMDs: {alt_spliced_count_with_tmds}", file=logfile)
    print(f"Number of TMDs predicted by TOPCONS2 in Ensembl CDS database: {num_unique_tmds}", file=logfile)
    print(f"Number of TMDs predicted by TOPCONS2 in filtered Ensembl data: {num_tmd_entries_in_filtered_data}", file=logfile)
    print(f"Separate count of transcripts in TOPCONS data: {num_transcripts_with_tmds}", file=logfile)
    print(f"Number of unique Ensembl gene IDs encountered in Ensembl CDS database: {num_unique_genes}", file=logfile)
    print(f"Number of unique Ensembl gene IDs with TMDs encountered in Ensembl CDS database: {num_unique_genes_with_tmds}", file=logfile)
    print(f"Number of unique HGNC symbols encountered in Ensembl CDS database: {num_unique_hgnc_symbols}", file=logfile)
    print(f"Number of unique HGNC symbols with TMDs encountered in Ensembl CDS database: {num_unique_hgnc_symbols_with_tmds}", file=logfile)
    print(f"Number of Harrington motifs identified: {num_harrington_motifs}", file=logfile)
    print(f"Number of Harrington motifs identified in canonical transcripts: {canonical_motifs}", file=logfile)
    print(f"Number of Harrington motifs identified in alternative isoform transcripts: {alt_motifs}", file=logfile)

    logfile.close()

    output_list = [total_canonical_count,
                   total_alt_spliced_count,
                   total_count,
                   canonical_count_with_tmds,
                   alt_spliced_count_with_tmds,
                   total_count_with_tmds,
                   num_unique_tmds,
                   num_tmd_entries_in_filtered_data,
                   num_transcripts_with_tmds,
                   num_unique_genes,
                   num_unique_genes_with_tmds,
                   num_unique_hgnc_symbols,
                   num_unique_hgnc_symbols_with_tmds]

    return output_list




def parse_tmd_csv(path_tmd_csv):
    tmd_data = {}
    with open(path_tmd_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            transcript_id = row['seq_name']
            tmd_data.setdefault(transcript_id, []).append({
                'adjusted_tmd_seq': row['adjusted_tmd_seq'],
                'len_adjusted_tmd': int(row['len_adjusted_tmd']),
                'adjusted_tmd_start_position': int(row['adjusted_tmd_start_position']),
                'adjusted_tmd_end_position': int(row['adjusted_tmd_end_position']),
                'adjusted_tmd_delta_G': float(row['adjusted_tmd_delta_G']),
                'adjusted_tmd_upstream_loop_length': int(row['adjusted_tmd_upstream_loop_length']),
                'adjusted_tmd_downstream_loop_length': int(row['adjusted_tmd_downstream_loop_length']),
                'tmd_type': row['tmd_type'],  # Add tmd_type to the data structure
                'tmd_orientation': row['tmhmm_predicted_tmd_topology'] # Add tmd orientation to data structure
            })
    return tmd_data


def get_slipsite_candidates(filtered_list, tmd_data, all_ss_candidates, candidate_filter=lambda s: True):
    """
    For each transcript in filtered_list, scans the CDS (split into codons) for candidate slip–sites 
    (i.e. 7-nt sequences obtained from a 3-codon window starting at the 3rd nucleotide of the first codon)
    that satisfy candidate_filter. For each candidate, it considers all upstream TMDs (from tmd_data) and computes:
      - gap_size: candidate codon index minus the TMD end (0-indexed)
      - deviation: gap_size minus 45 (i.e. how far off from the ideal 45 codon gap)
      - min_distance: same as deviation (for backward compatibility)
      - ss_coordinate_codon: candidate codon number (1-indexed)
      - ss_coordinate_nuc: candidate nucleotide coordinate (ss_coordinate_codon*3 + 2)
      - chr_coordinate: computed from transcript locus (transcript_data[3])
    Also returns the chosen TMD’s info.
    
    Returns a list of candidate dictionaries with keys:
      'transcript_id', 'canonical_status', 'slipsite_seq', 'slipsite_position',
      'gap_size', 'deviation', 'min_distance', 'ss_coordinate_codon', 'ss_coordinate_nuc',
      'chr_coordinate', 'tmd_info'
    """
    global friction_data, scoring_sys

    import re

    candidates = []
    for transcript_data in filtered_list:
        transcript_id = transcript_data[1]
        # Clean and uppercase the nucleotide sequence
        seq = re.sub("[^A-Za-z]", "", transcript_data[5]).upper()
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
        canonical_status = transcript_data[6]
        if transcript_id not in tmd_data:
            continue

        # Build list of upstream TMDs as tuples: (tm_end (0-indexed), tmd_info)
        tmd_list = tmd_data[transcript_id]
        upstream_tmds = []
        for tmd in tmd_list:
            tm_end = tmd['adjusted_tmd_end_position'] - 1
            upstream_tmds.append((tm_end, tmd))
        
        # For each possible candidate slip–site
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            slipsite_seq = heptamer_long[2:]  # last 7 nt of the window

            if slipsite_seq in friction_data:
                friction_score = friction_data[slipsite_seq]
            else:
                friction_score = frameshift.measure_slipsite_window(slipsite_seq, scoring_sys)

            if not candidate_filter(slipsite_seq):
                continue
            
            # For every upstream TMD (where tm_end < candidate index i), compute gap and deviation
            valid_candidates = []
            for (tm_end, tmd) in upstream_tmds:
                if tm_end < i:
                    gap_size = i - tm_end       # gap in codons
                    deviation = gap_size - 45    # how far off from ideal 45 codons
                    # We only consider candidates where deviation > -45 (i.e. gap_size > 0)
                    if deviation > -45:
                        valid_candidates.append((tm_end, tmd, gap_size, deviation))
            if not valid_candidates:
                continue
            
            # Choose the upstream TMD with the smallest absolute deviation from 45
            chosen_tm_end, chosen_tmd, chosen_gap_size, chosen_deviation = min(
                valid_candidates, key=lambda x: abs(x[3])
            )
            # For backward compatibility, set min_distance equal to the deviation.
            min_distance = chosen_deviation
            ss_coordinate_codon = i + 1         # 1-indexed codon number
            #ss_coordinate_nuc = ss_coordinate_codon * 3 + 2  # as in original code
            ss_coordinate_nuc = (i * 3) + 2 + 1 # 0-based gets to NNX N, add 2 gets to X, + 1 switches to 1-based indexing
            
            # Compute chromosomal coordinate from the transcript locus information.
            locus_coord_parts = transcript_data[3].split(':')
            try:
                locus_coord_start = int(locus_coord_parts[2])
            except ValueError:
                locus_coord_start = int(locus_coord_parts[3])
            chr_coordinate = locus_coord_start + ss_coordinate_nuc
            
            candidate = {
                'transcript_id': transcript_id,
                'canonical_status': canonical_status,
                'slipsite_seq': slipsite_seq,
                'slipsite_position': i,      # 0-indexed codon position
                'gap_size': chosen_gap_size,
                'deviation': chosen_deviation,
                'min_distance': min_distance,
                'ss_coordinate_codon': ss_coordinate_codon,
                'ss_coordinate_nuc': ss_coordinate_nuc,
                'chr_coordinate': chr_coordinate,
                'tmd_info': chosen_tmd,
                'friction_score': friction_score
            }
            candidates.append(candidate)

    if all_ss_candidates == 2:
        output_candidates = best_ss_per_tmd(candidates)
    else:
        output_candidates = candidates

    return output_candidates



def best_ss_per_tmd(candidates):
    # Now sweep through candidates and group the ones with a common transcript_id, TMD
    # Then find the downstream slip site closest to 45  
    output_candidates = []

    # First create a dictionary binning candidates in terms of their shared TMD
    tmd_grouped_candidates = {}
    for cand in candidates:
        transcript_id = cand['transcript_id']
        # Retrieve TMD info from candidate (nearest upstream TMD)
        tmd = cand['tmd_info']
        tm_start = tmd['adjusted_tmd_start_position']
        tm_end = tmd['adjusted_tmd_end_position']
        slipsite_seq = cand['slipsite_seq']
        gap_size = cand['gap_size']
        deviation = cand['deviation']
        tmd_attributes = [transcript_id, tm_start, tm_end]
        tmd_key = frozenset(tmd_attributes)
        if tmd_key not in tmd_grouped_candidates:
            tmd_grouped_candidates[tmd_key] = []
            tmd_grouped_candidates[tmd_key].append(cand)
        else:
            tmd_grouped_candidates[tmd_key].append(cand)
    
    # Next, iterate through candidate bins and identify the candidate with the smallest deviation. Add this to output
    for tmd_key, list_of_candidates in tmd_grouped_candidates.items():
        sorted_group = sorted(
            list_of_candidates,
            key=lambda x: (
                abs(x['deviation']),
                x['friction_score'],
                x['deviation']
            )
        )
        best_candidate = sorted_group[0]
        output_candidates.append(best_candidate)

    return output_candidates

import os, sqlite3, gzip, re, tempfile, gc, uuid
from pyteomics import parser

# Precompiled once
BOUNDARY_RX = re.compile(r'\bframe_transition=(\d+)\b')
FSTYPE_RX   = re.compile(r'\bfs_type=([^\|]+)\b', re.IGNORECASE)
TERM_RX     = re.compile(r'\bterminal_stop_codon=(yes|no)\b', re.IGNORECASE)

def _pick_tmp_dir(preferred: str = None) -> str:
    """Return a writable temp directory, preferring `preferred` if possible."""
    candidates = []
    if preferred:
        candidates.append(preferred)
    # Respect TMPDIR if set, then system temp, then cwd/_tmp
    env_tmp = os.environ.get('TMPDIR') or os.environ.get('TMP') or os.environ.get('TEMP')
    for p in (env_tmp, tempfile.gettempdir(), os.path.join(os.getcwd(), "_tmp")):
        if p:
            candidates.append(p)
    # pick the first we can create/write
    for d in candidates:
        try:
            os.makedirs(d, exist_ok=True)
            test_path = os.path.join(d, f".wtest_{uuid.uuid4().hex}")
            with open(test_path, "wb") as _:
                pass
            os.remove(test_path)
            return d
        except Exception:
            continue
    # last resort: current directory
    return os.getcwd()

def _sqlite_conn(db_path):
    # isolation_level=None → autocommit; we'll manage BEGIN/COMMIT ourselves
    conn = sqlite3.connect(db_path, isolation_level=None)
    cur = conn.cursor()
    # speed/robustness pragmas
    cur.execute("PRAGMA journal_mode=WAL;")
    cur.execute("PRAGMA synchronous=NORMAL;")
    cur.execute("PRAGMA temp_store=MEMORY;")
    cur.execute("PRAGMA mmap_size=134217728;")
    # schema
    cur.execute("""
        CREATE TABLE IF NOT EXISTS peptides (
            peptide TEXT PRIMARY KEY,
            header  TEXT
        );
    """)
    # If you also want a full mapping table in SQLite, add it—but you already stream map_tsv.
    return conn

def write_peptide_fasta_streaming(
    full_length_records,
    out_fa,
    map_tsv,
    digest_opts: DigestOptions,
    tmp_dir=None,
    commit_every=20000
):
    """
    Memory-safe, FragPipe-consistent digest writer.
    - Fully enzymatic only (semi_enzymatic comes from digest_opts)
    - Uses SQLite to dedupe peptides and store one representative header
    - Streams the full mapping to map_tsv(.gz)
    """
    rules = build_enzyme_rules(digest_opts)
    opts = digest_opts

    # ---- tmp sqlite path (writable) ----
    base_tmp = _pick_tmp_dir(tmp_dir)
    db_path = os.path.join(base_tmp, f"peptmp_dedupe_{os.getpid()}_{uuid.uuid4().hex}.sqlite")

    conn = _sqlite_conn(db_path)
    cur = conn.cursor()

    # Small in-RAM batch for representative headers
    # We INSERT OR IGNORE to keep the first header seen per peptide
    rep_batch = []  # list[(pep, hdr)]
    rep_stmt = "INSERT OR IGNORE INTO peptides (peptide, header) VALUES (?, ?);"

    n_rows = 0  # mapping lines written
    n_rep  = 0  # attempted rep inserts

    opener = gzip.open if map_tsv.endswith('.gz') else open
    with opener(map_tsv, 'wt') as map_out:
        for parent_hdr, prot_seq in full_length_records:
            parent_seq_U = prot_seq.upper()

            # NME variants
            for vidx, parent_seq in enumerate(nme_variants(parent_seq_U, opts.nme_mode)):
                parent_hdr_eff = parent_hdr if vidx == 0 else (parent_hdr + '|nme_variant=1')

                # header fields
                m = BOUNDARY_RX.search(parent_hdr)
                zero_prefix = int(m.group(1)) if m else len(parent_seq)

                mfs = FSTYPE_RX.search(parent_hdr)
                fs_type = (mfs.group(1).lower() if mfs else '0_frame')
                mt  = TERM_RX.search(parent_hdr)
                term_yes = (mt and mt.group(1).lower() == 'yes')

                is_shift_product = ('-1' in fs_type) or ('-2' in fs_type)
                ends_at_fs_stop  = term_yes and is_shift_product
                last_idx = len(parent_seq) - 1

                for enz, rule in rules.items():
                    # iterate peptides; regex=True to respect our regex rules
                    peps_iter = parser.cleave(
                        parent_seq,
                        rule,
                        missed_cleavages=opts.missed,
                        min_length=opts.min_len,
                        max_length=opts.max_len,
                        regex=True
                    )
                    for pep in peps_iter:
                        L = len(pep)
                        # (already length filtered by pyteomics, but keep the guard)
                        if L < opts.min_len or L > opts.max_len:
                            continue

                        # enumerate all occurrences to set novelty + coords (deterministic)
                        start = 0
                        while True:
                            idx = parent_seq.find(pep, start)
                            if idx == -1:
                                break
                            start = idx + 1
                            s = idx
                            e = idx + L - 1

                            if s < zero_prefix <= e:
                                novelty = 'transitional'
                            elif s >= zero_prefix:
                                novelty = 'frameshift'
                            else:
                                if e == last_idx and ends_at_fs_stop:
                                    if '-2' in fs_type:
                                        novelty = 'zero_frame_with_minus_two_stop'
                                    elif '-1' in fs_type:
                                        novelty = 'zero_frame_with_minus_one_stop'
                                    else:
                                        novelty = 'zero_frame'
                                else:
                                    novelty = 'zero_frame'

                            child_hdr = f'{parent_hdr_eff}|enz={enz}|digest={s}-{e}|novelty={novelty}'
                            pepU = pep.upper()

                            # (1) stream mapping line
                            map_out.write(f'{pepU}\t{child_hdr}\n')
                            n_rows += 1

                            # (2) buffer representative header insert
                            rep_batch.append((pepU, child_hdr))
                            n_rep += 1

                            # commit batch
                            if (n_rep % commit_every) == 0:
                                cur.execute("BEGIN;")
                                cur.executemany(rep_stmt, rep_batch)
                                cur.execute("COMMIT;")
                                rep_batch.clear()
                                gc.collect()

        # flush remaining buffered inserts
        if rep_batch:
            cur.execute("BEGIN;")
            cur.executemany(rep_stmt, rep_batch)
            cur.execute("COMMIT;")
            rep_batch.clear()

    # Emit peptide FASTA from SQLite
    with gzip.open(out_fa, 'wt') as fh:
        for pep, hdr in cur.execute("SELECT peptide, header FROM peptides;"):
            fh.write('>' + hdr + '\n')
            for i in range(0, len(pep), 60):
                fh.write(pep[i:i+60] + '\n')

    conn.close()
    try:
        os.remove(db_path)
    except OSError:
        pass

    print(f'Peptide FASTA written:  {out_fa}')
    print(f'Mapping table written: {map_tsv}')



def write_fasta(records, fh, width=60):
    """
    records : iterable of (header_without_gt, seq) pairs
    fh      : *text* file handle (plain or gzip opened in 'wt')
    """
    for hdr, seq in records:
        if hdr.startswith('>'):
            raise ValueError(f"Header already contains '>' : {hdr[:40]}…")
        fh.write('>' + hdr + '\n')
        for i in range(0, len(seq), width):
            fh.write(seq[i:i+width] + '\n')


from pyteomics import parser


def build_enzyme_rules(opts: DigestOptions):
    """
    Return dict[name -> regex rule] that matches FragPipe settings you provided.
    All rules are regexes; we’ll call parser.cleave(..., regex=True).
    """
    # Trypsin family
    trypsin_core = r'(?<=[KR])(?!P)'   # KR, no cut before P  (FragPipe “trypsin”)
    stricttryps  = r'(?<=[KR])'        # KR, no exception     (FragPipe “stricttrypsin”)

    # Chymotrypsin
    chymo_core = r'(?<=[FYW])(?!P)'    # FYW, not before P
    chymo_low  = r'(?<=[FYWL])(?!P)'   # FYWL, not before P
    chymo_rule = chymo_low if opts.chymo_specificity == "low" else chymo_core

    # Glu-C (FragPipe: Cuts DE, No cuts P -> not before P)
    gluc_rule = r'(?<=[DE])(?!P)' if opts.gluc_mode == "E_or_D" else r'(?<=[E])(?!P)'

    # Asp-N, Lys-C, Lys-N
    aspn  = r'(?=D)'                   # before D
    lysc  = r'(?<=K)(?!P)'             # after K, not before P (FragPipe showed No cuts P)
    lys_n = r'(?=K)'                   # before K (no exception in FragPipe)

    # Choose trypsin mode
    tryps_rule = stricttryps if opts.trypsin_mode == 'stricttrypsin' else trypsin_core

    return {
        'trypsin': tryps_rule,
        'lys-c'  : lysc,
        'gluc'   : gluc_rule,
        'aspn'   : aspn,
        'chymo'  : chymo_rule,
        'lys-n'  : lys_n,
    }

def nme_variants(seq: str, mode: str) -> list:
    """
    Return a list of parent protein sequences including NME variants as requested.
    mode: "off" | "bio" | "always"
    bio: excise M only if position 2 is small/uncharged (A,C,G,P,S,T,V)
    """
    seq = seq or ""
    if not seq or seq[0] != 'M':
        return [seq]
    if mode == "off":
        return [seq]
    if mode == "always":
        return [seq, seq[1:]]
    # biological heuristic
    if len(seq) >= 2 and seq[1] in {'A','C','G','P','S','T','V'}:
        return [seq, seq[1:]]
    return [seq]

def cleave_semi_enzymatic(seq, rule, missed, min_len, max_len):
    """
    Enumerate semi-enzymatic peptides by expanding around the canonical cut positions.
    """
    # canonical (fully-enzymatic) peptides → infer cut sites
    peps_full = list(parser.cleave(seq, rule,
                                   missed_cleavages=missed,
                                   min_length=None, max_length=None,
                                   regex=True))
    cut_sites = {0}
    end_sites = set()
    for p in peps_full:
        i = seq.find(p, 0)
        while i != -1:
            cut_sites.add(i)
            end_sites.add(i + len(p))
            i = seq.find(p, i+1)
    cut_sites.add(len(seq)); end_sites.add(len(seq))
    cut_sites = sorted(cut_sites); end_sites = sorted(end_sites)

    L = len(seq)
    emitted = set()
    # start at canonical start; end anywhere
    for s in cut_sites:
        for e in range(s + min_len, min(L, s + max_len) + 1):
            emitted.add(seq[s:e])
    # start anywhere; end at canonical end
    for e in end_sites:
        s_min = max(0, e - max_len)
        s_max = max(0, e - min_len)
        for s in range(s_min, s_max + 1):
            emitted.add(seq[s:e])
    return emitted


def write_peptide_fasta(full_length_records,
                        out_fa,
                        map_tsv,
                        digest_opts: DigestOptions):
    """
    All thresholds and flags come from digest_opts.
    """
    import gzip, re, collections
    opts = digest_opts
    rules = build_enzyme_rules(opts)

    boundary_rx = re.compile(r'\bframe_transition=(\d+)\b')
    fs_type_rx  = re.compile(r'\bfs_type=([^\|]+)\b', re.IGNORECASE)
    term_rx     = re.compile(r'\bterminal_stop_codon=(yes|no)\b', re.IGNORECASE)

    seq2repr = {}
    seq2all  = collections.defaultdict(list)

    for parent_hdr, prot_seq in full_length_records:
        for variant_idx, parent_seq in enumerate(nme_variants(prot_seq.upper(), opts.nme_mode)):
            parent_hdr_eff = parent_hdr if variant_idx == 0 else (parent_hdr + '|nme_variant=1')

            m = boundary_rx.search(parent_hdr)
            zero_prefix = int(m.group(1)) if m else len(parent_seq)

            mfs = fs_type_rx.search(parent_hdr)
            fs_type = (mfs.group(1).lower() if mfs else '0_frame')
            mt  = term_rx.search(parent_hdr)
            term_yes = (mt and mt.group(1).lower() == 'yes')

            is_shift_product = ('-1' in fs_type) or ('-2' in fs_type)
            ends_at_fs_stop  = term_yes and is_shift_product
            last_idx = len(parent_seq) - 1

            for enz, rule in rules.items():
                # fully-enzymatic
                full_set = set(parser.cleave(
                    parent_seq,
                    rule,
                    missed_cleavages=opts.missed,
                    min_length=None, max_length=None,
                    regex=True
                ))
                # choose peptide set
                if opts.semi_enzymatic:
                    peps = cleave_semi_enzymatic(parent_seq, rule, opts.missed, opts.min_len, opts.max_len)
                else:
                    peps = {p for p in full_set if opts.min_len <= len(p) <= opts.max_len}

                for pep in peps:
                    L = len(pep)
                    if not (opts.min_len <= L <= opts.max_len):
                        continue

                    # enumerate all occurrences to set novelty deterministically
                    start = 0
                    while True:
                        idx = parent_seq.find(pep, start)
                        if idx == -1:
                            break
                        start = idx + 1
                        s = idx
                        e = idx + L - 1

                        if s < zero_prefix <= e:
                            novelty = 'transitional'
                        elif s >= zero_prefix:
                            novelty = 'frameshift'
                        else:
                            if e == last_idx and ends_at_fs_stop:
                                if '-2' in fs_type:
                                    novelty = 'zero_frame_with_minus_two_stop'
                                elif '-1' in fs_type:
                                    novelty = 'zero_frame_with_minus_one_stop'
                                else:
                                    novelty = 'zero_frame'
                            else:
                                novelty = 'zero_frame'

                        child_hdr = f'{parent_hdr_eff}|enz={enz}|digest={s}-{e}|novelty={novelty}'
                        key = pep.upper()
                        if key not in seq2repr:
                            seq2repr[key] = child_hdr
                        seq2all[key].append(child_hdr)

    with gzip.open(out_fa, 'wt') as fh:
        for pep, hdr in seq2repr.items():
            fh.write('>' + hdr + '\n')
            for i in range(0, len(pep), 60):
                fh.write(pep[i:i+60] + '\n')

    with gzip.open(map_tsv, 'wt') as fh:
        for pep, hdrs in seq2all.items():
            for h in hdrs:
                fh.write(f'{pep}\t{h}\n')

    print(f'Peptide FASTA written:  {out_fa}')
    print(f'Mapping table written: {map_tsv}  '
          f'({len(seq2repr):,} unique peptides; {sum(map(len,seq2all.values())):,} total headers)')




def write_peptide_fasta_old(full_length_records,
                        out_fa,
                        map_tsv,
                        missed=2,
                        min_len=6, 
                        max_len=50):
    """
    full_length_records : list[(header, protein_seq)]  from find_fs_peptide_fastas
    out_fa  : gz-compressed peptide FASTA written here
    map_tsv : gzipped 2-col mapping  <peptide>\t<header>
      (one row per parent header)

    New novelty classes added:
      - zero_frame_with_minus_one_stop
      - zero_frame_with_minus_two_stop
    """

    from pyteomics import parser
    import gzip, re, collections

    # ----- enzyme rules (loose) -----
    t = parser.expasy_rules
    high = t['chymotrypsin high specificity']
    low  = t['chymotrypsin low specificity']
    tryps = t['trypsin']
    tryps_exc = t['trypsin_exception']
    gluc  = t['glutamyl endopeptidase']
    staph = t['staphylococcal peptidase i']
    aspn  = t['asp-n']
    lysc  = t['lysc']
    lys_n = r'(?=K)(?!P)'  # Lys-N

    rules = {
        'trypsin': f'({tryps})|({tryps_exc})',
        'lys-c'  : lysc,
        'gluc'   : f'({gluc})|({staph})',   # E or D
        'aspn'   : aspn,
        'chymo'  : f'({high})|({low})',     # high + low
        'lys-n'  : lys_n,
    }

    # ----- parent-header field parsers -----
    boundary_rx = re.compile(r'\bframe_transition=(\d+)\b')
    fs_type_rx  = re.compile(r'\bfs_type=([^\|]+)\b', re.IGNORECASE)
    term_rx     = re.compile(r'\bterminal_stop_codon=(yes|no)\b', re.IGNORECASE)

    # Dictionaries
    seq2repr = {}           # peptide -> representative header (str)
    seq2all = {}          # peptide -> list of *all* headers

    for parent_hdr, prot_seq in full_length_records:
        prot_seq = prot_seq.upper()
        # (1) locate 0→shift AA boundary (length of 0-frame prefix)
        m = boundary_rx.search(parent_hdr)
        zero_prefix = int(m.group(1)) if m else len(prot_seq)

        # (2) parent fs_type and whether end is a stop
        mfs = fs_type_rx.search(parent_hdr)
        fs_type = (mfs.group(1).lower() if mfs else '0_frame')
        mt  = term_rx.search(parent_hdr)
        term_yes = (mt and mt.group(1).lower() == 'yes')

        # convenience flags
        is_shift_product = ('-1' in fs_type) or ('-2' in fs_type)
        ends_at_fs_stop  = term_yes and is_shift_product
        last_idx = len(prot_seq) - 1  # 0-based index of final residue

        # (3) digest and classify peptides
        for enz, rule in rules.items():
            for pep in parser.cleave(prot_seq, rule, missed_cleavages=missed):
                if not (min_len <= len(pep) <= max_len):
                    continue

                # novelty class based on coordinates within parent
                start = prot_seq.find(pep)
                if start == -1:
                    # shouldn't happen with pyteomics fragments, but be defensive
                    continue
                end = start + len(pep) - 1

                if start < zero_prefix <= end:
                    novelty = 'transitional'
                elif start >= zero_prefix:
                    novelty = 'frameshift'
                else:
                    # wholly before the transition → usually zero_frame,
                    # but upgrade if peptide ends exactly at the protein end
                    # and that end is caused by a frameshift stop.
                    if end == last_idx and ends_at_fs_stop:
                        if '-2' in fs_type:
                            novelty = 'zero_frame_with_minus_two_stop'
                        elif '-1' in fs_type:
                            novelty = 'zero_frame_with_minus_one_stop'
                        else:
                            # should not occur, but fall back cleanly
                            novelty = 'zero_frame'
                    else:
                        novelty = 'zero_frame'

                # build child header (append our fields)
                child_hdr = (
                    f'{parent_hdr}'
                    f'|enz={enz}'
                    f'|digest={start}-{end}'
                    f'|novelty={novelty}'
                )

                key = pep.upper()
                if key not in seq2repr:
                    seq2repr[key] = child_hdr
                    seq2all[key]  = [child_hdr]
                else:
                    seq2all[key].append(child_hdr)

    # Write peptide FASTA of uniques
    with gzip.open(out_fa, 'wt') as fh:
        for pep, hdr in seq2repr.items():
            fh.write('>' + hdr + '\n')
            for i in range(0, len(pep), 60):
                fh.write(pep[i:i+60] + '\n')

    # Write peptide→parent mapping (one row per parent header)
    with gzip.open(map_tsv, 'wt') as fh:
        for pep, hdrs in seq2all.items():
            for h in hdrs:
                fh.write(f'{pep}\t{h}\n')

    print(f'Peptide FASTA written:  {out_fa}')
    print(f'Mapping table written: {map_tsv}  '
          f'({len(seq2repr):,} unique peptides; {sum(map(len,seq2all.values())):,} total headers)')


from collections import defaultdict
from textwrap import wrap
import gzip, re, os.path as op

from collections import defaultdict, Counter
import gzip, re

def _zopen(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def write_fasta(records, fh, width=60):
    for hdr, seq in records:
        if hdr.startswith('>'):
            raise ValueError(f"Header already contains '>' : {hdr[:40]}…")
        fh.write('>' + hdr + '\n')
        for i in range(0, len(seq), width):
            fh.write(seq[i:i+width] + '\n')


def build_unique_and_fsonly(src_fa, unique_fa, fs_only_fa, mapping_tsv):
    """
    1) unique_fa  – deduplicated peptide FASTA (one representative header per peptide)
    2) fs_only_fa – peptides that are *never* zero_frame across all parents
                    (i.e., only frameshift-evidence classes below)

    Frameshift-evidence classes:
        - 'frameshift'
        - 'transitional'
        - 'zero_frame_with_minus_one_stop'
        - 'zero_frame_with_minus_two_stop'
    """
    hdr_re = re.compile(r'novelty=([^|\n]+)')     # capture novelty tag

    pep2parents = defaultdict(list)               # peptide -> [headers…]
    pep2novelty = defaultdict(set)                # peptide -> {novelty classes…}

    # ingest
    with _zopen(src_fa, 'rt') as fh:
        header, seq_chunks = None, []
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    seq = ''.join(seq_chunks).strip().upper()
                    pep2parents[seq].append(header.lstrip('>').rstrip())
                    m = hdr_re.search(header)
                    if m:
                        pep2novelty[seq].add(m.group(1).lower())
                header, seq_chunks = line.rstrip('\n'), []
            else:
                seq_chunks.append(line.strip())
        # flush last
        if header is not None:
            seq = ''.join(seq_chunks).strip().upper()
            pep2parents[seq].append(header.lstrip('>').rstrip())
            m = hdr_re.search(header)
            if m:
                pep2novelty[seq].add(m.group(1).lower())

    # evidence sets
    EVIDENCE_CLASSES = {
        'frameshift',
        'transitional',
        'zero_frame_with_minus_one_stop',
        'zero_frame_with_minus_two_stop',
    }

    # Write TSV mapping and prep FASTAs
    uniq_records, fs_only_records = [], []
    evidence_any_ct = 0

    with gzip.open(mapping_tsv, 'wt') as mapfh:
        mapfh.write('peptide\tparents\n')
        for pep, parents in pep2parents.items():
            master = f"{parents[0]}|sources={len(parents)}"  # no leading '>'
            uniq_records.append((master, pep))
            mapfh.write(f"{pep}\t" + ';'.join(parents) + '\n')

            novset = pep2novelty.get(pep, set())
            # broader "evidence_any": peptide has any evidence class, regardless of also being zero_frame somewhere
            if novset & EVIDENCE_CLASSES:
                evidence_any_ct += 1
            # strict fs-only FASTA: peptide is never zero_frame anywhere
            if 'zero_frame' not in novset:
                fs_only_records.append((master, pep))

    # Write FASTAs
    with gzip.open(unique_fa, 'wt') as fh:
        write_fasta(uniq_records, fh)
    with gzip.open(fs_only_fa, 'wt') as fh:
        write_fasta(fs_only_records, fh)

    # Logs
    total_unique = len(pep2parents)
    fs_only_ct   = len(fs_only_records)

    print(f"{total_unique:,} unique peptide sequences written → {unique_fa}")
    print(f"Strict frameshift-only subset (never zero_frame): {fs_only_ct:,} → {fs_only_fa}")
    print(f"Frameshift-evidence peptides (any evidence class, even if also zero_frame elsewhere): {evidence_any_ct:,}")
    print(f"Mapping written → {mapping_tsv}")


'''
def translate_from(nt: str, start_nt: int) -> str:
    gc_dict = genetic_code_func()
    aa = []
    for i in range(start_nt, len(nt)-2, 3):
        codon = nt[i:i+3]
        aa.append(gc_dict.get(codon, 'X'))
        if aa[-1] == '*':
            break
    return ''.join(aa)
'''



STOP_AA = '*'

def _clean_seq(nt: str) -> str:
    return ''.join(b for b in nt.upper() if b in 'ACGT')

def _translate_from(nt: str, start_nt: int, gc_dict) -> Tuple[str, bool]:
    """Translate starting at a 0-based nucleotide offset (step 3) until first stop.
       Returns (aa_without_stop, terminated_bool)."""
    aa = []
    terminated = False
    # guard for out-of-range starts
    if start_nt < 0:
        start_nt = 0
    for i in range(start_nt, len(nt) - 2, 3):
        codon = nt[i:i+3]
        aa_sym = gc_dict.get(codon, 'X')
        aa.append(aa_sym)
        if aa_sym == STOP_AA:
            aa.pop()      # do not include the stop in the returned sequence
            terminated = True
            break
    return ''.join(aa), terminated

def _zero_prefix_through(nt: str, codon_index_inclusive_0based: int, gc_dict) -> str:
    """Translate 0-frame from AUG through and including the given 0-based codon index,
       without stopping early on internal stops (we need the prefix length in AA)."""
    end_nt_excl = (codon_index_inclusive_0based + 1) * 3
    aa = []
    for i in range(0, min(end_nt_excl, len(nt) - 2), 3):
        aa.append(gc_dict.get(nt[i:i+3], 'X'))
    return ''.join(aa)


def compute_dual_m1_fs_tail(nt_seq, ss_codon_1based, gc_dict=None):
    """
    Compute the dual −1 (two-tRNA) frameshift peptide *tail only*:
      - 0-frame up through YYZ (not returned here)
      - then start translation at (end of YYZ) − 1 nt (−1 frame)
      - translate until stop or transcript end

    Returns:
        fs_tail_with_optional_star, terminated ("yes" or "no")
    """
    if gc_dict is None:
        gc_dict = genetic_code_func()

    nt = _clean_seq(nt_seq)

    # codon/window geometry — keep identical to find_fs_peptide_fastas
    win0 = ss_codon_1based - 1           # 0-based codon index of NNX
    i_yyz = win0 + 2                     # codon index of YYZ
    nt_end_YYZ = (i_yyz + 1) * 3         # nucleotide index just after YYZ

    start_dual_m1 = nt_end_YYZ - 1       # one nt before next 0-frame codon → −1 frame

    # not enough nucleotides to even form one codon
    if start_dual_m1 > len(nt) - 3:
        return "", "no"

    # Look at the very first −1-frame codon
    first_codon = nt[start_dual_m1:start_dual_m1+3]
    first_aa = gc_dict.get(first_codon, 'X')

    # If the very first −1 codon is a stop, represent the FS product as "*"
    if first_aa == STOP_AA:
        return "*", "yes"

    # Otherwise, let the generic translator handle the rest
    tail_dual_m1, term_dual_m1 = _translate_from(nt, start_dual_m1, gc_dict)

    # _translate_from already trimmed any internal stop; attach '*' if we saw one
    fs_peptide = tail_dual_m1 + ('*' if term_dual_m1 else '')
    terminated = 'yes' if term_dual_m1 else 'no'

    return fs_peptide, terminated



def find_fs_peptide_fastas(harrington_motifs_with_uniprot, fasta_path=None):
    """
    Input: list of motif rows with schema you showed (indices 0..31 used here).
    Output: list of (header:str, seq:str) for 4 products/motif:
            0_frame, dual_-1, hungry_-1, hungry_-2
    """

    fasta_full_length = []
    gc_dict = genetic_code_func()

    for entry in harrington_motifs_with_uniprot:
        # names & meta
        transcript_id   = entry[0]
        gene_id         = entry[1]
        gene_name       = entry[2]
        tmd_start       = entry[3]
        tmd_end         = entry[4]
        tmd_seq         = entry[6]
        ss_codon_1based = entry[9]                # first codon of NNX|XXY|YYZ window
        slip_heptamer   = entry[11]
        gap_size        = entry[12]
        nuc_seq_raw     = entry[19]               # full CDS (DNA)
        ensembl_uniprot = entry[31]
        # optional / may be absent depending on your run
        blast_uniprot   = entry[32] if len(entry) > 32 else 'NA'
        canonical_status= entry[26]
        tmd_dg          = entry[5]
        friction_score  = entry[21]
        stem_loop       = entry[23]

        # clean nucleotide sequence
        nt = _clean_seq(nuc_seq_raw)

        # codon/window geometry
        win0 = ss_codon_1based - 1            # 0-based index of NNX (window start)

        i_xxy = win0 + 1                      # 0-based index of XXY
        i_yyz = win0 + 2                      # 0-based index of YYZ

        nt_end_XXY = (i_xxy + 1) * 3          # end of XXY (exclusive)
        nt_end_YYZ = (i_yyz + 1) * 3          # end of YYZ (exclusive)
        
        
        # 0-frame prefixes (AA counts needed for headers)
        prefix_thru_XXY = _zero_prefix_through(nt, i_xxy, gc_dict)  # includes XXY
        prefix_thru_YYZ = _zero_prefix_through(nt, i_yyz, gc_dict)  # includes YYZ

        # verify the 9mer window matches the recorded heptamer
        # sanity check that this candidate’s 7-mer matches the CDS-derived 7-mer
        header_mismatch_tag = False
        if slip_heptamer:  # e.g., your entry[11]
            nine_mer = nt[win0*3 : win0*3 + 9]
            if len(nine_mer) == 9:
                hept_from_nt = nine_mer[2:]  # last 7 nt of this window
                if hept_from_nt != slip_heptamer:
                    # Don’t shift to ±1 here: adjacent windows are intentional candidates.
                    # Record a soft warning for debugging, but continue with this window.
                    # (Optional) add a tag to the header:
                    header_mismatch_tag = f"heptamer_mismatch={slip_heptamer}->{hept_from_nt}"
            else:
                header_mismatch_tag = "heptamer_mismatch=short_window"
        else:
            header_mismatch_tag = "heptamer_mismatch=none"
        
        # Comment this out when not testing
        if header_mismatch_tag != False:
            print(f"{header_mismatch_tag}")

        # --- build four products ---

        products = {}

        # (1) 0-frame full-length from AUG to first stop
        zero_full, zero_term = _translate_from(nt, 0, gc_dict)
        products['0_frame'] = (zero_full, zero_term, 'no_frameshift')

        # (2) dual −1 (two-tRNA): resume at (end YYZ) − 1 nt, keep 0-frame through YYZ
        start_dual_m1 = nt_end_YYZ - 1
        tail_dual_m1, term_dual_m1 = _translate_from(nt, start_dual_m1, gc_dict)
        products['dual_-1'] = (
            prefix_thru_YYZ + tail_dual_m1,
            term_dual_m1,
            str(len(prefix_thru_YYZ))
        )

        # (3) hungry −1: resume at (end XXY) − 1 nt, keep 0-frame through XXY
        start_hungry_m1 = nt_end_XXY - 1
        tail_hungry_m1, term_hungry_m1 = _translate_from(nt, start_hungry_m1, gc_dict)
        products['hungry_-1'] = (
            prefix_thru_XXY + tail_hungry_m1,
            term_hungry_m1,
            str(len(prefix_thru_XXY))
        )

        # (4) hungry −2: resume at (end XXY) − 2 nt, keep 0-frame through XXY
        seq_hungry_m2, term_hungry_m2 = '', False
        start_hungry_m2 = nt_end_XXY - 2
        if start_hungry_m2 >= 0:
            tail_hungry_m2, term_hungry_m2 = _translate_from(nt, start_hungry_m2, gc_dict)
            seq_hungry_m2 = prefix_thru_XXY + tail_hungry_m2
        products['hungry_-2'] = (
            seq_hungry_m2,
            term_hungry_m2,
            str(len(prefix_thru_XXY))
        )

        # --- headers & output (keeps your field names/order) ---
        # Note: you previously computed ss_coordinate_nuc as (ss_codon*3)+2 (0-based + 2) in discovery code;
        # if you need that exact value in the header, recompute here the same way:
        #ss_coordinate_nuc = ss_codon_1based * 3 + 2
        #ss_coordinate_nuc = entry[8] # use stored value rather than re-calculate # nucleotide coordinate of third X in 'NNX XXY YYZ'
        ss_coordinate_nuc = (win0 * 3) + 2 + 1 # nucleotide coordinate of first X in 'NNX XXY YYZ '

        base_header_fields = [
            f'ensembl_transcript_id={transcript_id}',
            f'ensembl_gene_id={gene_id}',
            f'gene_name={gene_name}',
            f'tmd_start={tmd_start}',
            f'tmd_end={tmd_end}',
            f'tmd_seq={tmd_seq}',
            f'ss_nuc_start={ss_coordinate_nuc}',
            f'ss_seq={slip_heptamer}',
            f'tmd-SS_gap_size={gap_size}',
            f'ensembl_uniprot_annotation={ensembl_uniprot}',
            f'blastp_uniprot_annotation={blast_uniprot}',
            f'transcript_canonical_status={canonical_status}',
            f'TMD_dG_app={tmd_dg:.2f}',
            f'SS_dG_fs={friction_score:.2f}',
            f'stem_loop_presence={stem_loop}'
        ]

        for kind in ('0_frame', 'dual_-1', 'hungry_-1', 'hungry_-2'):
            aa_seq, terminated, ftpos = products[kind]
            header_suffix = f'frame_transition={ftpos}'
            term_suffix   = 'terminal_stop_codon=yes' if terminated else 'terminal_stop_codon=no'
            header = '|'.join(base_header_fields) + f'|fs_type={kind}|{header_suffix}|{term_suffix}'

            # IMPORTANT: do not include '*' in the output AA (we already trimmed at stop)
            fasta_full_length.append((header, aa_seq))

    return fasta_full_length



def motif_search(filtered_list, tmd_data, slipsite_list, friction_data, gap_boundary_near, gap_boundary_far, scoring_sys, all_ss_candidates):
    """
    canonical_status = filtered_list[6]
    ensembl_uniprot_annotation = filtered_list[7]
    ensembl_transcript_type = filtered_list[8]
    Searches for Harrington motifs by first scanning each transcript for slip–site candidates
    using get_slipsite_candidates (with candidate_filter ensuring the slipsite is in slipsite_list).
    Then, for each candidate, enriches it with additional features to form a candidate entry with fields:
      [0] transcript_id
      [1] gene_id
      [2] gene_name
      [3] tm_start (1-indexed)
      [4] tm_end (1-indexed)
      [5] tm_dg
      [6] tm_domain_str (the protein segment corresponding to the TMD)
      [7] tm_codons_str (the coding sequence for the TMD)
      [8] ss_coordinate_nuc
      [9] ss_coordinate_codon
      [10] chr_coordinate
      [11] candidate_slipsite (the 7-nt slipsite)
      [12] gap_size (in codons)
      [13] deviation from ideal (gap_size - 45)
      [14] gap_nuc (nucleotide sequence between TMD and slip site)
      [15] gap_gc_content (% GC in gap)
      [16] gap_avg_codon_freq (average codon frequency in gap)
      [17] fs_peptide (frameshifted peptide downstream)
      [18] terminated (whether translation terminated with a stop)
      [19] transcript_seq (nucleotide sequence)
      [20] transcript_residues_str (protein sequence)
      [21] friction_score
      [22] slipsite_strict ("Strict" if slipsite in set_24_canonical_slipsites else "Loose")
      [23] stem_loop_found ("Yes" or "No")
      [24] pseudoknot_found ("Yes" or "No")
      [25] secondary_structure (predicted structure)
      [26] canonical_status_ensembl_annotation
      [27] tmd_type
      [28] tmd_orientation
      [29] transcript_type_ensembl_annotation
      [30] ensembl_transcript_aliases
      [31] uniprot_annotation_ensembl
    Also counts how many motifs fall within [gap_boundary_near, gap_boundary_far] and tallies trial counts.
    
    Returns:
        motif_list, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative
    """
    # Define strict slip sites
    set_24_canonical_slipsites = {'AAAAAAA','AAAAAAT','AAAAAAC','AAATTTA','AAATTTT','AAATTTC',
                                  'TTTAAAA','TTTAAAT','TTTAAAC','TTTTTTA','TTTTTTT','TTTTTTC',
                                  'CCCAAAA','CCCAAAT','CCCAAAC','CCCTTTA','CCCTTTT','CCCTTTC',
                                  'GGGAAAA','GGGAAAT','GGGAAAC','GGGTTTA','GGGTTTT','GGGTTTC'}
    
    # Get candidates that pass candidate_filter (only if in slipsite_list)
    candidates = get_slipsite_candidates(filtered_list, tmd_data, all_ss_candidates, candidate_filter=lambda s: s in slipsite_list)
    
    # Build a lookup for transcript info, codons, and protein sequence.
    transcript_info = {}
    transcript_codons_dict = {}
    transcript_protein_dict = {}
    for transcript_data in filtered_list:
        transcript_id = transcript_data[1]
        transcript_info[transcript_id] = {
            'gene_id': transcript_data[0],
            'gene_name': transcript_data[2],
            'locus': transcript_data[3],
            'seq': re.sub("[^a-zA-Z]", "", transcript_data[5]),
            'canonical_status': transcript_data[6],
            'ensembl_uniprot_annotation' : transcript_data[7],
            'ensembl_transcript_type': transcript_data[8],
            'ensembl_transcript_aliases': transcript_data[9]
        }
        seq_codons = __import__("textwrap").wrap(transcript_info[transcript_id]['seq'], 3)
        transcript_codons_dict[transcript_id] = seq_codons
        transcript_protein_dict[transcript_id] = ribosome(seq_codons, genetic_code_func())
    
    motif_list = []
    motifs_within_gap = 0
    num_trials_all = 0
    num_trials_canonical = 0
    num_trials_alternative = 0

    # Process each candidate
    for cand in candidates:
        transcript_id = cand['transcript_id']
        info = transcript_info[transcript_id]
        codons = transcript_codons_dict[transcript_id]
        residues = transcript_protein_dict[transcript_id]
        
        # Retrieve TMD info from candidate (nearest upstream TMD)
        tmd = cand['tmd_info']
        tm_start = tmd['adjusted_tmd_start_position']
        tm_end = tmd['adjusted_tmd_end_position']
        tm_dg = tmd['adjusted_tmd_delta_G']
        tmd_type = tmd.get('tmd_type', 'NA')
        # Get TMD protein segment and coding sequence:
        tm_protein = residues[tm_start - 1: tm_end]
        tm_domain_str = ''.join(tm_protein)
        tm_codons_str = ''.join(codons[tm_start - 1: tm_end])
        
        # Coordinates from candidate
        ss_coordinate_nuc = cand['ss_coordinate_nuc']
        ss_coordinate_codon = cand['ss_coordinate_codon']
        chr_coordinate = cand['chr_coordinate']
        candidate_slipsite = cand['slipsite_seq']
        gap_size = cand['gap_size']
        deviation = cand['deviation']
        
        # Determine gap nucleotide sequence from TMD end to candidate slipsite (1-indexed: use indices tm_end to ss_coordinate_codon-1)
        gap_codons = codons[tm_end: ss_coordinate_codon - 1] if ss_coordinate_codon - 1 > tm_end else []
        gap_nuc_raw = ''.join(gap_codons)
        gap_nuc = gap_nuc_raw[:-1] if gap_nuc_raw else ''
        gap_gc_content = get_GC(gap_nuc) if gap_nuc else 0
        gap_avg_codon_freq = calc_avg_freq(gap_codons) if gap_codons else 0
        
        # Count trial if transcript is long enough (we mimic original: if len(codons) >= tm_end + 45 + 2)
        if len(codons) >= tm_end + 45 + 2:
            num_trials_all += 1
            if info['canonical_status'] == 'Canonical':
                num_trials_canonical += 1
            else:
                num_trials_alternative += 1
        
        # Count motifs within user–defined gap boundaries (apply to gap_size)
        if gap_boundary_near <= gap_size <= gap_boundary_far:
            motifs_within_gap += 1
        
        # Compute friction score using the 9-mer window starting at candidate slipsite:
        candidate_window = ''.join(codons[cand['slipsite_position'] : cand['slipsite_position'] + 3])
        if candidate_slipsite in friction_data:
            friction_score = friction_data[candidate_slipsite]
        else:
            friction_score = frameshift.measure_slipsite_window(candidate_window, scoring_sys)
        
        # Compute downstream sequence for frameshift peptide
        # Compute -1 frameshifted segment of dual t-RNA frameshifting only;
        # for more products see find_fs_peptide_fastas function
        fs_peptide, terminated = compute_dual_m1_fs_tail(
            info['seq'],                 # full CDS (nucleotides)
            cand['ss_coordinate_codon']  # 1-based codon index of NNX (window start)
        )

        # Predict mRNA secondary structure in a downstream window (60 nt window starting at (ss_coordinate_codon*3)+5)
        downstream_window_start = (ss_coordinate_codon * 3) + 5
        downstream_window_end = min(len(info['seq']), downstream_window_start + 60)
        downstream_window = info['seq'][downstream_window_start:downstream_window_end]
        structure, mfe = predict_secondary_structure(downstream_window)
        stem_loop_found = "Yes" if detect_stem_loop(structure) else "No"
        pseudoknot_found = "Yes" if detect_pseudoknot(structure) else "No"
        
        # Determine slipsite strictness
        slipsite_strict = "Strict" if candidate_slipsite in set_24_canonical_slipsites else "Loose"
        
        # Build the candidate entry with all fields (following original order):
        candidate_entry = [
            transcript_id,                        # 0
            info['gene_id'],                      # 1
            info['gene_name'],                    # 2
            tm_start,                             # 3
            tm_end,                               # 4
            tm_dg,                                # 5
            tm_domain_str,                        # 6
            tm_codons_str,                        # 7
            ss_coordinate_nuc,                    # 8
            ss_coordinate_codon,                  # 9
            chr_coordinate,                       # 10
            candidate_slipsite,                   # 11
            gap_size,                             # 12
            deviation,                            # 13
            gap_nuc,                              # 14
            gap_gc_content,                       # 15
            gap_avg_codon_freq,                   # 16
            fs_peptide,                           # 17
            terminated,                           # 18
            info['seq'],                          # 19 transcript nucleotide sequence
            ''.join(residues),                    # 20 transcript protein sequence
            friction_score,                       # 21
            slipsite_strict,                      # 22
            stem_loop_found,                      # 23
            pseudoknot_found,                     # 24
            structure,                            # 25
            info['canonical_status'],             # 26
            tmd.get('tmd_type', 'NA'),            # 27
            tmd.get('tmd_orientation', 'NA'),            # 28
            info['ensembl_transcript_type'],      # 29
            info['ensembl_transcript_aliases'],   # 30
            info['ensembl_uniprot_annotation']    # 31
        ]
        motif_list.append(candidate_entry)
    

    #motif_list_all = motif_list

    motif_list_within_gap = []
    for entry in motif_list:
        gap_size = int(entry[12])
        if gap_size in range(gap_boundary_near, gap_boundary_far + 1):
            motif_list_within_gap.append(entry)

    return motif_list, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative, motif_list_within_gap



# Print motif results to CSV
def print_motif_csv(output_array, outname):
    with open(outname + '.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Updated headers to match the structure of candidate_entry and include Uniprot information
        csv_writer.writerow([
            'transcript_id', 'gene_id', 'gene_name', 'tm_start', 'tm_end', 'tm_dg', 'tm_domain', 'tm_codons',
            'ss_coord_nuc', 'ss_coord_codon', 'chr_coordinate', 'slipsite_seq', 'gap_size', 'location_ideality',
            'gap_nuc', 'gap_%gc', 'gap_avg_codon_freq', 'fs_peptide', 'terminated', 'transcript_seq_nuc',
            'transcript_seq_res', 'friction_score', 'slipsite_canonical?', 'stem_loop_found', 'pseudoknot_found', 'secondary_structure',
            'canonical_status', 'tmd_type', 'tmd_orientation', 'ensembl_transcript_type', 'ensembl_transcript_id_aliases', 'ensembl_uniprot_id', 'blasp_uniprot_id', 'blasp_similarity_score'
        ])
        
        # Write motif entries to the CSV
        for motif in output_array:
            csv_writer.writerow(motif)



def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()


def print_ensembl_uniprot_dict(ensembl_uniprot_dict, outname):
    ensembl_uniprot_list = [['ensembl_transcript_id', 'uniprot_accession_id', 'similarity_score']]
    for key, value in ensembl_uniprot_dict.items():
        ensembl_transcript_id = key
        [uniprot_accession_id, similarity_score] = value
        output_entry = [ensembl_transcript_id, uniprot_accession_id, similarity_score]
        ensembl_uniprot_list.append(output_entry)
    output_csv_name = outname+'_ensembl_uniprot_dict_full_blastp'
    print_results(ensembl_uniprot_list, output_csv_name)



def split_up_harrington_motifs(harrington_motifs, location_ideality_threshold, stem_loop_required=False):
    '''
    canonical_status = filtered_list[6]
    ensembl_uniprot_annotation = filtered_list[7]
    ensembl_transcript_type = filtered_list[8]
    Searches for Harrington motifs by first scanning each transcript for slip–site candidates
    using get_slipsite_candidates (with candidate_filter ensuring the slipsite is in slipsite_list).
    Then, for each candidate, enriches it with additional features to form a candidate entry with fields:
      [0] transcript_id
      [1] gene_id
      [2] gene_name
      [3] tm_start (1-indexed)
      [4] tm_end (1-indexed)
      [5] tm_dg
      [6] tm_domain_str (the protein segment corresponding to the TMD)
      [7] tm_codons_str (the coding sequence for the TMD)
      [8] ss_coordinate_nuc
      [9] ss_coordinate_codon
      [10] chr_coordinate
      [11] candidate_slipsite (the 7-nt slipsite)
      [12] gap_size (in codons)
      [13] deviation from ideal (gap_size - 45)
      [14] gap_nuc (nucleotide sequence between TMD and slip site)
      [15] gap_gc_content (% GC in gap)
      [16] gap_avg_codon_freq (average codon frequency in gap)
      [17] fs_peptide (frameshifted peptide downstream)
      [18] terminated (whether translation terminated with a stop)
      [19] transcript_seq (nucleotide sequence)
      [20] transcript_residues_str (protein sequence)
      [21] friction_score
      [22] slipsite_strict ("Strict" if slipsite in set_24_canonical_slipsites else "Loose")
      [23] stem_loop_found ("Yes" or "No")
      [24] pseudoknot_found ("Yes" or "No")
      [25] secondary_structure (predicted structure)
      [26] canonical_status_ensembl_annotation
      [27] tmd_type
      [28] tmd_orientation
      [29] transcript_type_ensembl_annotation
      [30] ensembl_transcript_aliases
      [31] uniprot_annotation_ensembl
      [32] uniprot_id_blasp   
      [33] blasp_similarity_score
    '''

    # Initialize strict subsets
    strict_all_motifs = []
    strict_canonical_motifs = []
    strict_alternative_motifs = []

    # Initialize loose subsets
    loose_all_motifs = []
    loose_canonical_motifs = []
    loose_alternative_motifs = []

    for entry in harrington_motifs:
        location_ideality = abs(int(entry[13]))
        slipsite_strict = entry[22]
        canonical_status = entry[26]
        terminated_status = entry[18]
        stem_loop_presence = entry[23]

        subset_dimensions = [slipsite_strict, canonical_status]

        if all([location_ideality <= location_ideality_threshold,
                stem_loop_required == True,
                stem_loop_presence == 'yes']) or all([stem_loop_required == False,
                                                      location_ideality <= location_ideality_threshold]):

            if subset_dimensions == ['Strict', 'Canonical']:
                strict_all_motifs.append(entry)
                strict_canonical_motifs.append(entry)
                loose_all_motifs.append(entry)
                loose_canonical_motifs.append(entry)
        
            if subset_dimensions == ['Strict', 'Alternative']:
                strict_all_motifs.append(entry)
                strict_alternative_motifs.append(entry)
                loose_all_motifs.append(entry)
                loose_alternative_motifs.append(entry)
        
            if subset_dimensions == ['Loose', 'Canonical']:
                loose_all_motifs.append(entry)
                loose_canonical_motifs.append(entry)
        
            if subset_dimensions == ['Loose', 'Alternative']:
                loose_all_motifs.append(entry)
                loose_alternative_motifs.append(entry)


    motif_subsets = [strict_all_motifs, strict_canonical_motifs, strict_alternative_motifs, loose_all_motifs, loose_canonical_motifs, loose_alternative_motifs]

    return motif_subsets


def configparser_function(config_file="config.ini"):
    """
    Loads configuration from an ini file using configparser.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Paths
    path_fasta = config.get('Paths', 'path_fasta')
    path_uniprot = config.get('Paths', 'path_uniprot')
    blastp_results_path = config.get('Paths', 'blastp_results_path')
    path_tmd_csv = config.get('Paths', 'path_tmd_csv')
    path_slipsites = config.get('Paths', 'path_slipsites')
    path_friction_csv = config.get('Paths', 'path_friction_csv')
    path_ensembl_features = config.get('Paths', 'path_ensembl_features')
    
    # Parameters
    gap_near = config.getint('Parameters', 'gap_near')
    gap_far = config.getint('Parameters', 'gap_far')
    outname = config.get('Parameters', 'outname')
    scoring_sys = config.getint('Parameters', 'scoring_sys')  # Will return 1 or 2
    all_ss_candidates = config.getint('Parameters', 'all_ss_candidates') # Will return 1 or 2
    
    # Return all values in the same order as interactive_configuration
    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys, all_ss_candidates


# Main configuration function
def interactive_configuration():
    path_fasta = input("Please enter a path to an Ensembl-format CDS database: ")
    path_uniprot = input("Please enter a path to an Excel sheet containing Uniprot sequence data: ")
    blastp_results_path = input("Please enter a path to the output of the blastp search matching Ensembl transcripts to Uniprot sequences: ")
    path_tmd_csv = input("Please enter a path to the von Heijne optimized TMD data CSV: ")
    path_slipsites = input("Please enter a path to a list of slip sites to search for: ")
    path_friction_csv = input("Please enter a path to the friction scores CSV: ")
    path_ensembl_features = input("Please enter a path to the TSV containing features (HGNC, Uniprot, etc) from Ensembl: ")
    gap_near = int(input("Please enter a lower bound on the distance range in the sequence upstream from the slippery sequence to be searched (in codons): "))
    gap_far = int(input("Please enter an upper bound on the distance range in the sequence upstream from the slippery sequence to be searched (in codons): "))
    outname = input("Please enter a name for the output: ")
    while True:
        print("Which scoring system do you want to use")
        print("[1] Use number of hydrogen-bonds formed and broken (simple scoring)")
        print("[2] Use Turner nearest neighbor parameters")
        scoring_sys = int(input("Please enter [1] or [2]: "))
        if scoring_sys in (1, 2):
            break
        else:
            print("There has been an error. Please re-enter.")

    while True:
        print("Do you want all slippery sequences downstream of a TMD as unique candidates?")
        print("[1] For every TMD, create a candidate for every slippery sequence within the window.")
        print("[2] For every TMD, create a candidate with the slippery sequence closest to the ideal spacing.")
        all_ss_candidates = int(input("Please enter [1] or [2]: "))
        if all_ss_candidates in (1, 2):
            break
        else:
            print("There has been an error. Please re-enter.")
    

    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys, all_ss_candidates

def write_reference_proteomes(filtered_list, outname):

    canonical_cds = [entry for entry in filtered_list if entry[6]=='Canonical']


    write_reference_proteome(canonical_cds, outname+"_Homo_sapiens.GRCh38.filtered.pep.canonical.fa.gz")
    write_reference_proteome(filtered_list, outname+"_Homo_sapiens.GRCh38.filtered.pep.all.fa.gz")


def write_reference_proteome(transcripts, out_fasta):
    # This is the indexing of the output
    # output_line = [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len(sequence_str), sequence_str, canonical_status, uniprot_id, transcript_type, alias_transcript_ids]
    #                    0                 1                2           3                 4                 5             6               7               8               9


    with _zopen(out_fasta, 'wt') as fh:
        for gene_id, tx_id, gene_name, locus, length, nuc_seq, canonical_status, uniprot_id, *_ in transcripts:
            # translate
            codons = [nuc_seq[i:i+3] for i in range(0, len(nuc_seq)-2, 3)]
            prot = ''.join(ribosome(codons, genetic_code_func()))
            # drop terminal '*' if present
            if prot.endswith('*'):
                prot = prot[:-1]
            # header—pick whatever you like here
            header = f">transcript_id={tx_id}|gene_id={gene_id}|{gene_name}|canonical_status={canonical_status}|uniprot_match={uniprot_id}"
            fh.write(header + "\n")
            for line in wrap(prot, 60):
                fh.write(line + "\n")

def is_canonical(header: str) -> bool:

    return "transcript_canonical_status=Canonical" in header


def write_fs_chimeric_fastas(all_entries, outname):
    """
    Build two gzipped FASTAs:
      - out_all_fa: every chimera product you generate
      - out_canon_fa: only those from transcripts with canonical_status == 'Canonical'
    """
    # split
    all_chimeras = all_entries
    canon_chimeras = [(h,s) for (h,s) in all_entries if is_canonical(h)]

    out_all_fa = f'{outname}_all_full_length_chimeras.fa.gz'
    out_canon_fa = f'{outname}_canonical_full_length_chimeras.fa.gz'

    import gzip
    if out_all_fa:
        with gzip.open(out_all_fa, 'wt') as fh:
            write_fasta(all_chimeras, fh)
        print(f"All chimeras written tp {out_all_fa}")

    if out_canon_fa:
        with gzip.open(out_canon_fa, 'wt') as fh:
            write_fasta(canon_chimeras, fh)
        print(f"Canonical-only chimeras written to {out_canon_fa}")
    
    return canon_chimeras

def main():
    
    global friction_data, scoring_sys

    if len(sys.argv) < 2:
        print("Usage: python script.py <config_file.ini>")
        config_mode_flag = int(input("Do you ant to revert to interactive mode? [1] yes or [2] no"))
        if config_mode_flag == 1:
            path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys, all_ss_candidates = interactive_configuration()
        else:
            sys.exit(1)
    
    config_file = sys.argv[1]

    path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys, all_ss_candidates = configparser_function(config_file)

    fasta_list = parse_fasta(path_fasta)

    stop_filtered_fasta_list = stop_codon_filter(fasta_list)
    slipsite_list = import_text_file(path_slipsites)
    tmd_data = parse_tmd_csv(path_tmd_csv)
    print("Size of tmd_data: ", len(tmd_data))
    friction_data = parse_friction_csv(path_friction_csv)

    # Create dictionary with ensembl features
    ensembl_features_dict = parse_ensembl_uniprot_canon(path_ensembl_features)

    print("Processing transcripts and matching to Uniprot sequences...")
    filtered_list, canonical_count, alt_spliced_count = fasta_transcript_filter(stop_filtered_fasta_list, ensembl_features_dict)

    import pickle
    filtered_list_pickle = f'{outname}_filtered_list.pkl'
    with open(filtered_list_pickle, 'wb') as f:
        pickle.dump(filtered_list, f)
    print(f"filtered_list successfully picked to {filtered_list_pickle}")

    print("Printing reference proteomes from Ensembl filtered sequences.")
    write_reference_proteomes(filtered_list, outname)


    ensembl_uniprot_dict = process_transcripts_with_uniprot(filtered_list, blastp_results_path, tmd_data, outname)
    print_ensembl_uniprot_dict(ensembl_uniprot_dict, outname)

    # Perform the Harrington motif search
    print("Search for Harrington motifs...")
    all_motifs, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative, harrington_motifs = motif_search(filtered_list, tmd_data, slipsite_list, friction_data, gap_near, gap_far, scoring_sys, all_ss_candidates)
    
    print("Associaing Harrington motifs with Uniprot sequences")
    harrington_motifs_with_uniprot = append_uniprot_accession_to_motif(harrington_motifs, ensembl_uniprot_dict)
    all_motifs_with_uniprot = append_uniprot_accession_to_motif(all_motifs, ensembl_uniprot_dict)


    # Generatign fasta with full-length chimeric sequence (0 reading frame to -1 reading frame)
    print("Building chimeric peptides that transition to alternative reading frame for each motif and writing to a fasta")

    fasta_full_length_all = find_fs_peptide_fastas(harrington_motifs_with_uniprot)

    fasta_full_length_canonical = write_fs_chimeric_fastas(fasta_full_length_all, outname)

 
    print("Digesting chimeras …")

    dig_opts = DigestOptions(
        gluc_mode="E_or_D",       # Glu-C DE, !P
        chymo_specificity="low",  # FLWY, !P
        trypsin_mode="stricttrypsin",  # or "trypsin" if that’s what they ran
        semi_enzymatic=False,
        missed=2,
        min_len=7,
        max_len=50,
        nme_mode="always"
    )

    write_peptide_fasta_streaming(
        fasta_full_length_all,
        out_fa   = f'{outname}_digested_fs_pep_all.fa.gz',
        map_tsv  = f'{outname}_peptide2parent_all.tsv.gz',
        digest_opts = dig_opts,
        commit_every = 20000
    )

    print("Building unique and frameshift-only subsets …")
    build_unique_and_fsonly(
        src_fa      = f'{outname}_digested_fs_pep_all.fa.gz',
        unique_fa   = f'{outname}_digest_unique_all.fa.gz',
        fs_only_fa  = f'{outname}_digest_nonzeroframe_only_all.fa.gz',
        mapping_tsv = f'{outname}_peptide_source_mapping_all.tsv.gz')

    # CANONICAL-only chimeras
    write_peptide_fasta_streaming(
        fasta_full_length_canonical,
        out_fa   = f'{outname}_digested_fs_pep_canonical.fa.gz',
        map_tsv  = f'{outname}_peptide2parent_canonical.tsv.gz',
        digest_opts = dig_opts,
        commit_every = 20000
    )

    print("Building unique and frameshift-only subsets …")
    build_unique_and_fsonly(
        src_fa      = f'{outname}_digested_fs_pep_canonical.fa.gz',
        unique_fa   = f'{outname}_digest_unique_canonical.fa.gz',
        fs_only_fa  = f'{outname}_digest_nonzeroframe_only_canonical.fa.gz',
        mapping_tsv = f'{outname}_peptide_source_mapping_canonical.tsv.gz')

    print('Done; peptide FASTA ready for MS search.')


    print("Obtaining some basic information about transcript and motif counts")
    harrington_output_list = get_canonical_count(filtered_list, tmd_data, harrington_motifs, outname+'_harrington_motifs')
    all_output_list = get_canonical_count(filtered_list, tmd_data, all_motifs, outname+'_all_motifs')

    print_motif_csv(harrington_motifs_with_uniprot, outname+'_harrington_motifs')
    print_motif_csv(all_motifs_with_uniprot, outname+'_all_motifs')
    
    print(f"Basic statistics for Harrington motifs (candidates within gap boundaries): ")
    print(f"Total canonical transcripts: {canonical_count}")
    print(f"Total alternatively-spliced transcripts: {alt_spliced_count}")
    print(f"Number of motifs within the gap range ({gap_near}-{gap_far} codons): {motifs_within_gap}")



if __name__ == '__main__':

    main()
