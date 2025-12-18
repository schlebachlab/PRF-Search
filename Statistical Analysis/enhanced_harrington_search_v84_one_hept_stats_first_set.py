#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to search transcriptome databases (Ensembl CDS) for Harrington motifs
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
from collections import Counter
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
                'tmd_type': row['tmd_type']  # Add tmd_type to the data structure
            })
    return tmd_data


def get_slipsite_candidates(filtered_list, tmd_data, candidate_filter=lambda s: True):
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
            ss_coordinate_nuc = ss_coordinate_codon * 3 + 2  # as in original code
            
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


    #return candidates # For V74
    return output_candidates # For V73


def motif_search(filtered_list, tmd_data, slipsite_list, friction_data, gap_boundary_near, gap_boundary_far, scoring_sys):
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
      [28] transcript_type_ensembl_annotation
      [29] ensembl_transcript_aliases
      [30] uniprot_annotation_ensembl
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
    candidates = get_slipsite_candidates(filtered_list, tmd_data, candidate_filter=lambda s: s in slipsite_list)
    
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
        downstream_seq = ''.join(codons[cand['slipsite_position']:])
        if len(downstream_seq) >= 3:
            downstream_codons_fs = __import__("textwrap").wrap(downstream_seq[2:-1], 3)
            downstream_fs_peptide_full = ribosome(downstream_codons_fs, genetic_code_func())
            fs_peptide_list = []
            for residue in downstream_fs_peptide_full:
                if residue != '*':
                    fs_peptide_list.append(residue)
                else:
                    fs_peptide_list.append(residue)
                    break
            fs_peptide = ''.join(fs_peptide_list)
            terminated = 'yes' if fs_peptide and fs_peptide[-1] == '*' else 'no'
        else:
            fs_peptide = ''
            terminated = 'no'

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
            info['ensembl_transcript_type'],      # 28
            info['ensembl_transcript_aliases'],   # 29
            info['ensembl_uniprot_annotation']    # 30
        ]
        motif_list.append(candidate_entry)
    

    #motif_list_all = motif_list

    motif_list_within_gap = []
    for entry in motif_list:
        gap_size = int(entry[12])
        if gap_size in range(gap_boundary_near, gap_boundary_far + 1):
            motif_list_within_gap.append(entry)

    return motif_list, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative, motif_list_within_gap


def analyze_slipsite_to_tmd_distances(filtered_list, tmd_data, slipsite_list, outname):
    strict_slip_sites = {'AAAAAAA', 'AAAAAAT', 'AAAAAAC', 'AAATTTA', 'AAATTTT', 'AAATTTC',
                         'TTTAAAA', 'TTTAAAT', 'TTTAAAC', 'TTTTTTA', 'TTTTTTT', 'TTTTTTC',
                         'CCCAAAA', 'CCCAAAT', 'CCCAAAC', 'CCCTTTA', 'CCCTTTT', 'CCCTTTC',
                         'GGGAAAA', 'GGGAAAT', 'GGGAAAC', 'GGGTTTA', 'GGGTTTT', 'GGGTTTC'}
    # Accept candidates if they are in either the loose set or the strict set.
    candidates = get_slipsite_candidates(filtered_list, tmd_data,
                                          candidate_filter=lambda s: s in slipsite_list or s in strict_slip_sites)
    
    strict_all, strict_canonical, strict_alternative = [], [], []
    loose_all, loose_canonical, loose_alternative = [], [], []
    for cand in candidates:
        if cand['slipsite_seq'] in strict_slip_sites:
            strict_all.append(cand['min_distance'])
            if cand['canonical_status'] == 'Canonical':
                strict_canonical.append(cand['min_distance'])
            else:
                strict_alternative.append(cand['min_distance'])
        else:
            loose_all.append(cand['min_distance'])
            if cand['canonical_status'] == 'Canonical':
                loose_canonical.append(cand['min_distance'])
            else:
                loose_alternative.append(cand['min_distance'])
    
    distance_data = {
        "strict_all": strict_all,
        "strict_canonical": strict_canonical,
        "strict_alternative": strict_alternative,
        "loose_all": loose_all,
        "loose_canonical": loose_canonical,
        "loose_alternative": loose_alternative
    }
    
    # Write CSV file with the categorized distances
    with open(f"{outname}_v55_slipsite_tmd_distances_strict_and_loose_slipsite_distribution.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Distance_from_TMD', 'Category'])
        for d in strict_all:
            writer.writerow([d, 'strict_all'])
        for d in strict_canonical:
            writer.writerow([d, 'strict_canonical'])
        for d in strict_alternative:
            writer.writerow([d, 'strict_alternative'])
        for d in loose_all:
            writer.writerow([d, 'loose_all'])
        for d in loose_canonical:
            writer.writerow([d, 'loose_canonical'])
        for d in loose_alternative:
            writer.writerow([d, 'loose_alternative'])
    
    return distance_data



def analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list):
    # Only keep candidates that match the provided slipsite_list.
    candidates = get_slipsite_candidates(filtered_list, tmd_data, candidate_filter=lambda s: s in slipsite_list)
    distances = [cand['min_distance'] for cand in candidates]
    num_ideal = sum(1 for cand in candidates if cand['min_distance'] == 0)
    return distances, num_ideal


def analyze_slipsite_to_tmd_distances_contingency_table(filtered_list, tmd_data, slipsite_list, all_heptamers):
    counts = {'slippery_ideal': 0, 'non_slippery_ideal': 0, 'slippery_other': 0, 'non_slippery_other': 0}
    
    for transcript_data in filtered_list:
        transcript_id = transcript_data[1]
        sequence_raw = transcript_data[5]
        sequence = ''.join(filter(str.isalpha, sequence_raw)).upper()
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        if transcript_id not in tmd_data:
            continue
        tmd_positions = [tmd['adjusted_tmd_end_position'] - 1 for tmd in tmd_data[transcript_id]]
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            heptamer = heptamer_long[2:]
            is_slippery = heptamer in slipsite_list
            upstream = [tm for tm in tmd_positions if tm < i]
            if not upstream:
                continue
            distances = [(i - tm) - 45 for tm in upstream]
            valid = [d for d in distances if d > -45]
            if not valid:
                continue
            min_distance = min(valid, key=lambda d: abs(d))
            if min_distance == 0:
                if is_slippery:
                    counts['slippery_ideal'] += 1
                else:
                    counts['non_slippery_ideal'] += 1
            else:
                if is_slippery:
                    counts['slippery_other'] += 1
                else:
                    counts['non_slippery_other'] += 1
    return counts


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
            'canonical_status', 'tmd_type', 'ensembl_transcript_type', 'ensembl_transcript_id_aliases', 'ensembl_uniprot_id', 'blasp_uniprot_id', 'blasp_similarity_score'
        ])
        
        # Write motif entries to the CSV
        for motif in output_array:
            csv_writer.writerow(motif)





def draw_random_heptamers(num_slipsite_samples, shuffled_list_heptamers):

    random_heptamers = set()
    
    while len(random_heptamers) < num_slipsite_samples:
        random_heptamer_candidate = random.choice(shuffled_list_heptamers)
        random_heptamers.add(random_heptamer_candidate)
    
    return list(random_heptamers)


def import_uniprot_fasta(path_fasta):
    fasta_generator = load_fasta_data(path_fasta)
    raw_fasta_output = []
    for header, sequence in fasta_generator:
        sequence_str = ''.join(sequence)
        output_line = [header, sequence_str]
        raw_fasta_output.append(output_line)

    uniprot_fasta_data = []
    for entry in raw_fasta_output:
        header = entry[0]
        uniprot_seq = entry[1]
        uniprot_id = header.split('|')[1]
        gene_name = header.split('=')[-3].split()[0]
        output_line = [uniprot_id, gene_name, uniprot_seq]
        uniprot_fasta_data.append(output_line)
    return uniprot_fasta_data


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



def load_uniprot_excel_to_dict(excel_path, key_column='Entry'):
    """
    Loads a Uniprot Excel file and converts it into a dictionary with the specified key column.
    
    Parameters:
    excel_path (str): Path to the Excel file containing Uniprot data.
    key_column (str): The column to use as the key for the dictionary (e.g., 'Sequence', 'Entry').
    
    Returns:
    dict: A dictionary where keys are from the key_column and values are dictionaries of the other fields, 
          including the sequence itself.
    """
    # Load the Excel file into a Pandas DataFrame
    df = pd.read_excel(excel_path)

    # Fill any NaN values with empty strings to avoid issues with missing data
    df = df.fillna('')

    # Ensure the key_column exists in the DataFrame
    if key_column not in df.columns:
        raise ValueError(f"Key column '{key_column}' does not exist in the Excel file.")
    
    # Add the sequence as a value in the dictionary (even though it's the key)
    # Create the dictionary and make sure 'Sequence' is part of the values as well.
    uniprot_dict = df.set_index(key_column).T.to_dict('dict')

    # Explicitly add 'Sequence' into the dictionary values
    #for seq in uniprot_dict.keys():
    #    uniprot_dict[seq]['Sequence'] = seq  # Add the sequence as a value

    return uniprot_dict



def filter_non_canonical_transcripts(filtered_list):
    filtered_list_canonical_only = []
    filtered_list_alternative_only = []
    for entry in filtered_list:
        canonical_status = entry[6]
        if canonical_status == 'Canonical':
            filtered_list_canonical_only.append(entry)
        if canonical_status == 'Alternative':
            filtered_list_alternative_only.append(entry)
    return filtered_list_canonical_only, filtered_list_alternative_only


def extract_heptamers_for_prf(sequence_str):
    """
    Extract heptamers that could be involved in -1 PRF from a nucleotide sequence in codon space.
    The heptamer starts at the third base of the first codon and spans 7 nucleotides.

    Args:
        sequence_str: Nucleotide sequence (CDS) as a string.
        
    Returns:
        List of heptamers that could be involved in -1 PRF.
    """
    # Step 1: Split the nucleotide sequence into codons
    seq_codons = textwrap.wrap(sequence_str, 3)
    
    seq_heptamers = []
    
    # Step 2: Build the heptamers by starting at the 3rd base of the first codon
    for i in range(len(seq_codons) - 2):  # Need at least 3 codons to form a heptamer
        codon1 = seq_codons[i][2:]    # Take the 3rd base of the first codon
        codon2 = seq_codons[i+1]      # Take all of the second codon
        codon3 = seq_codons[i+2]      # Take all of the third codon
        
        # Combine them to form a heptamer (7 nucleotides)
        heptamer = codon1 + codon2 + codon3
        
        seq_heptamers.append(heptamer)
    
    return seq_heptamers




def calculate_heptamer_frequencies_prf(filtered_list, tmd_data):
    """
    Calculate the frequency of heptamers that could be involved in -1 PRF in the CDS sequences.

    Args:
        filtered_list: List of CDS sequences from the transcriptome.
        
    Returns:
        heptamer_frequencies: A dictionary of heptamer frequencies.
    """
    heptamer_counts = Counter()
    total_heptamers = 0
    
    list_of_heptamers = []

    filtered_list_with_tmds = []
    for entry in filtered_list:
        transcript_id = entry[1]
        if transcript_id in tmd_data:
            filtered_list_with_tmds.append(entry)

    # Iterate through each transcript's CDS sequence
    for transcript_data in filtered_list_with_tmds:
        sequence_str = transcript_data[5]  # CDS sequence (nucleotide space)
        
        # Extract heptamers for -1 PRF
        seq_heptamers = extract_heptamers_for_prf(sequence_str)
        for heptamer in seq_heptamers:
            list_of_heptamers.append(heptamer)
        
        # Count the heptamers
        for heptamer in seq_heptamers:
            heptamer_counts[heptamer] += 1
            total_heptamers += 1
    
    set_unique_heptamers = list(set(list_of_heptamers))
    shuffled_list_of_heptamers = list_of_heptamers.copy()
    random.shuffle(shuffled_list_of_heptamers)

    # Convert counts to frequencies
    heptamer_frequencies = {heptamer: count / total_heptamers for heptamer, count in heptamer_counts.items()}
    
    return heptamer_frequencies, set_unique_heptamers, shuffled_list_of_heptamers




def calculate_p_null_slippery(heptamer_frequencies, slippery_sequences):
    """
    Calculate the probability of finding a slippery sequence by summing their frequencies.

    Args:
        heptamer_frequencies: A dictionary of heptamer frequencies.
        slippery_sequences: A set of slippery sequences to consider.
        
    Returns:
        p_null: The probability of observing a slippery sequence at random.
    """
    p_null = sum(heptamer_frequencies.get(seq, 0) for seq in slippery_sequences)
    return p_null




def process_iteration_for_tmd_shuffling(args):
    filtered_list, tmd_data, slipsite_list, seed_i = args

    random.seed(seed_i)
    np.random.seed(seed_i)

    shuffled_tmd_data = shuffle_tmd_positions(filtered_list, tmd_data)
    shuffled_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, shuffled_tmd_data, slipsite_list)
    return shuffled_distances, num_ideal_distances



def run_tmd_shuffling_test_parallel(filtered_list, tmd_data, slipsite_list, base_seed, num_iterations=16384, num_processes=32):
    """
    Perform the TMD shuffling test in parallel and compare distances to the original TMD distances.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary of TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of shuffling iterations to run (default is 1000).
        num_processes: Number of parallel processes (default is 4).
    
    Returns:
        shuffled_distances: A list of distances between shuffled TMDs and slip sites (one list for each iteration).
    """
    
    # Prepare the arguments for each process
    args = [(filtered_list, tmd_data, slipsite_list, base_seed + i) for i in range(num_iterations)]

    # Use multiprocessing Pool to parallelize the iterations
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_iteration_for_tmd_shuffling, args) 

    # Unpack the results
    shuffled_distances_list = []
    num_ideal_distances_list = []

    for shuffled_distances, num_ideal_distances in results:
        shuffled_distances_list.append(shuffled_distances)
        num_ideal_distances_list.append(num_ideal_distances)

    return shuffled_distances_list, num_ideal_distances_list






def shuffle_tmd_positions(filtered_list, tmd_data):
    """
    Shuffle the positions of TMDs within each transcript while maintaining the constraints:
    - Preserve the number of TMDs.
    - Preserve the length of each TMD.
    - Preserve the order of TMDs.
    - Ensure TMDs do not overlap.
    - Maintain at least minimal_spacing codons between TMDs.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        minimal_spacing: Minimum number of codons required between TMDs (default is 3).
        
    Returns:
        shuffled_tmd_data: A new dictionary of shuffled TMD positions for each transcript.
    """
    # SET MINIMAL SPACING HERE
    minimal_spacing = 3
    # ADJUST MINIMAL SPACING IF DESIRED

    shuffled_tmd_data = {}
    
    # Build a dictionary with transcript sequences (for convenience)
    transcript_sequences = {entry[1]: entry[5] for entry in filtered_list}
    
    # Iterate over each transcript
    for transcript_id, tmd_list in tmd_data.items():
        if transcript_id not in transcript_sequences:
            continue
        
        sequence = transcript_sequences[transcript_id]
        seq_length = len(sequence) // 3  # Convert nucleotide sequence length to codon length
        num_tmds = len(tmd_list)
        
        # Get the lengths of each TMD
        tmd_lengths = [
            tmd['adjusted_tmd_end_position'] - tmd['adjusted_tmd_start_position'] + 1
            for tmd in tmd_list
        ]
        total_tmd_length = sum(tmd_lengths)
        
        # Calculate the minimum total gap required between TMDs
        total_min_gap = (num_tmds - 1) * minimal_spacing
        total_gap_space = seq_length - total_tmd_length
        
        # Check if it's possible to place TMDs under constraints
        extra_space = total_gap_space - total_min_gap
        if extra_space < 0:
            # Cannot fit TMDs under constraints
            continue
        
        # Distribute extra_space randomly among the gaps
        num_gaps = num_tmds + 1  # Gaps before, between, and after TMDs
        minimal_gaps = [0] + [minimal_spacing] * (num_tmds - 1) + [0]
        
        # Randomly distribute extra_space among the gaps
        delta_gaps = [0] * num_gaps
        if extra_space > 0:
            random_points = sorted(random.sample(range(extra_space + num_gaps - 1), num_gaps - 1))
            delta_gaps = [random_points[0]] + \
                         [random_points[i] - random_points[i - 1] - 1 for i in range(1, num_gaps - 1)] + \
                         [extra_space + num_gaps - 2 - random_points[-1]]
        
        # Calculate the actual gaps
        gaps = [minimal_gaps[i] + delta_gaps[i] for i in range(num_gaps)]
        
        # Determine the start positions for each TMD
        start_positions = []
        current_position = gaps[0]
        for i in range(num_tmds):
            tmd_start = current_position
            tmd_end = tmd_start + tmd_lengths[i] - 1
            start_positions.append((tmd_start, tmd_end, tmd_list[i]['tmd_type']))
            current_position = tmd_end + 1 + gaps[i + 1]
        
        # Build the new TMD list with shuffled positions
        new_tmd_list = []
        for tmd_info, (new_start, new_end, tmd_type) in zip(tmd_list, start_positions):
            new_tmd_list.append({
                'adjusted_tmd_start_position': new_start,
                'adjusted_tmd_end_position': new_end,
                'tmd_type': tmd_type,
            })
        
        # Add the shuffled TMD positions to the new dictionary
        shuffled_tmd_data[transcript_id] = new_tmd_list
    
    return shuffled_tmd_data


def iqr(arr):
    q75, q25 = np.percentile(arr, [75, 25])
    return q75 - q25

def perm_test(real_val, null_vals):
    mu = null_vals.mean()
    sigma = null_vals.std(ddof=1)
    z = (real_val - mu) / sigma
    pz = 2 * stats.norm.sf(abs(z))
    extreme = np.sum(np.abs(null_vals - mu) >= abs(real_val - mu))
    pe = (extreme + 1) / (len(null_vals) + 1)
    return z, pz, pe, mu, sigma, extreme




def compare_tmd_shuffle_to_observed(filtered_list, tmd_data, slipsite_list, base_seed, num_iterations=16384):
    """
    Compare real slipsite-to-TMD distances to shuffled TMD distances and compute Z-scores or p-values.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of shuffling iterations to run (default is 1000).
    
    Returns:
        z_score: Z-score comparing real distances to shuffled distances.
        p_value: Two-sided p-value indicating whether real distances are significantly different from shuffled.
    """
    
    # Step 1: Calculate real distances between slip sites and TMDs
    print("Calculating real slipsite-to-TMD distances...")
    real_distances, num_ideal_distances_real = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list)
    
    # Step 2: Run the TMD shuffling test to generate shuffled distances
    print(f"Running TMD shuffling test with {num_iterations} iterations...")
    shuffled_distances_list, num_ideal_distances_shuffled = run_tmd_shuffling_test_parallel(filtered_list, tmd_data, slipsite_list, base_seed, num_iterations=16384, num_processes=32)

    print("Comparing real distances to shuffled distances...")

    # Compute the real distances mean
    real_mean_distance = np.mean(real_distances)

    # Compute the null means (one per iteration)
    null_means = np.array([np.mean(distances) for distances in shuffled_distances_list])


    # Obtain z_perm, p_z, p_perm, mu_null, and sigma_null for permutation test
    z_score_permutation, p_z, p_perm, mu_null, sigma_null, extreme_count = perm_test(real_mean_distance, null_means)
    print(f'**********\n')
    print(f'PERMUTATION TEST RESULTS')    

    print(f'Unshuffled mean distance: {real_mean_distance}')
    print(f'Shuffled mean distance (mu_null): {mu_null}')
    print(f'Shuffled mean distance sigma (sigma_null): {sigma_null}\n')
    
    print(f'Z score: {z_score_permutation}')
    print(f'p_z: {p_z}\n')
    
    print(f'Extreme count: {extreme_count}')
    print(f'p_perm: {p_perm}\n')


    # Step 3: Compare real distances to shuffled distances using Z-score and p-value

    (ks_test_stat, ks_p_value, real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std, mw_stat, mw_p_value, real_zero_count, real_total, p_null_binom, binom_p_value,
        real_range_count, p_null_range_binom, binom_p_value_range,
        a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
        a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range) = compare_distances(real_distances, shuffled_distances_list)
    
    return (z_score_permutation, 
            p_z, 
            real_mean_distance, mu_null, sigma_null, 
            extreme_count, p_perm, 
            ks_test_stat, ks_p_value, 
            real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std,
            mw_stat, mw_p_value, 
            real_zero_count, real_total, p_null_binom, binom_p_value,
            real_range_count, p_null_range_binom, binom_p_value_range,
            num_ideal_distances_real, num_ideal_distances_shuffled, real_distances, shuffled_distances_list,
            a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
            a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range)



## Random heptamer distance test


def process_iteration_for_heptamer_sampling(args):
    filtered_list, tmd_data, slipsite_list, shuffled_list_of_heptamers, seed_i = args
    num_slipsite_samples = len(slipsite_list)
    
    # Set a new seed for this iteration
    random.seed(seed_i)
    np.random.seed(seed_i)

    # Randomly sample heptamers in proportion to their frequency
    random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers)
    
    # Calculate distances for the randomly sampled heptamers
    sampled_heptamer_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, random_heptamers)
    
    return sampled_heptamer_distances, num_ideal_distances


def run_heptamer_sampling_test_parallel(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, base_seed, num_iterations=16384, num_processes=32):
    """
    Perform a heptamer sampling test in parallel and compare their distances to TMDs.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of sampling iterations to run (default is 1000).
        num_processes: Number of parallel processes (default is 4).
    
    Returns:
        sampled_heptamer_distances: A list of distances between sampled heptamers and TMDs (one list for each iteration).
    """
    
    # Prepare the arguments for each process
    args = [(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, base_seed + i) for i in range(num_iterations)]

    # Use multiprocessing Pool to parallelize the iterations
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_iteration_for_heptamer_sampling, args)

    # Unpack the results
    sampled_heptamer_distances_list = []
    num_ideal_distances_list = []

    for sampled_heptamer_distances, num_ideal_distances in results:
        sampled_heptamer_distances_list.append(sampled_heptamer_distances)
        num_ideal_distances_list.append(num_ideal_distances)


    return sampled_heptamer_distances_list, num_ideal_distances_list


def compare_heptamer_sampling_to_observed(filtered_list, tmd_data, slipsite_list, base_seed, shuffled_heptamer_list, num_iterations=16384):
    """
    Compare real slipsite-to-TMD distances to distances from random heptamer sampling.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of sampling iterations to run (default is 1000).
    
    Returns:
        z_score: Z-score comparing real distances to sampled distances.
        p_value: P-value indicating whether the real distances are significantly different from sampled distances.
    """
    
    # Step 1: Calculate real distances between slip sites and TMDs
    print("Calculating real slipsite-to-TMD distances...")
    real_distances, num_ideal_distances_real = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list)
    
    # Step 2: Run the heptamer sampling test to generate random distances
    print(f"Running heptamer sampling test with {num_iterations} iterations...")
    sampled_heptamer_distances_list, num_ideal_distances_sampled = run_heptamer_sampling_test_parallel(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, base_seed, num_iterations=16384, num_processes=32)

    print("Comparing real distances to shuffled distances...")

    # Compute the real distances mean
    real_mean_distance = np.mean(real_distances)

    # Compute the null means (one per iteration)
    null_means = np.array([np.mean(distances) for distances in sampled_heptamer_distances_list])


    # Obtain z_perm, p_z, p_perm, mu_null, and sigma_null for permutation test
    z_score_permutation, p_z, p_perm, mu_null, sigma_null, extreme_count = perm_test(real_mean_distance, null_means)
    print(f'**********\n')
    print(f'PERMUTATION TEST RESULTS')    

    print(f'Unshuffled mean distance: {real_mean_distance}')
    print(f'Shuffled mean distance (mu_null): {mu_null}')
    print(f'Shuffled mean distance sigma (sigma_null): {sigma_null}\n')
    
    print(f'Z score: {z_score_permutation}')
    print(f'p_z: {p_z}\n')
    
    print(f'Extreme count: {extreme_count}')
    print(f'p_perm: {p_perm}\n')


    # Step 3: Compare real distances to shuffled distances using Z-score and p-value

    (ks_test_stat, ks_p_value, real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std, mw_stat, mw_p_value, real_zero_count, real_total, p_null_binom, binom_p_value,
        real_range_count, p_null_range_binom, binom_p_value_range,
        a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
        a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range) = compare_distances(real_distances, sampled_heptamer_distances_list)
    
    return (z_score_permutation, 
            p_z, 
            real_mean_distance, mu_null, sigma_null, 
            extreme_count, p_perm, 
            ks_test_stat, ks_p_value, 
            real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std,
            mw_stat, mw_p_value, 
            real_zero_count, real_total, p_null_binom, binom_p_value,
            real_range_count, p_null_range_binom, binom_p_value_range,
            num_ideal_distances_real, num_ideal_distances_sampled, real_distances, sampled_heptamer_distances_list,
            a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
            a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range)


import numpy as np
from scipy.stats import norm, ks_2samp, mannwhitneyu, binom_test

def get_distribution_stats(data):
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    std = np.std(data)
    return median, iqr, std

def mann_whitney_restricted(real, shuffled, lower_bound=-10, upper_bound=10):
    real_restricted = [d for d in real if lower_bound <= d <= upper_bound]
    shuffled_restricted = [d for d in shuffled if lower_bound <= d <= upper_bound]
    stat, p_val = mannwhitneyu(real_restricted, shuffled_restricted, alternative='two-sided')
    return stat, p_val

def binomial_test_at_zero(real, shuffled, tolerance=0):
    real_zero_count = sum(1 for d in real if abs(d) <= tolerance)
    shuffled_zero_count = sum(1 for d in shuffled if abs(d) <= tolerance)
    total_shuffled = len(shuffled)
    p_null = shuffled_zero_count / total_shuffled if total_shuffled > 0 else 0.0
    real_total = len(real)
    # Using binom_test (or binomtest in newer versions)
    p_value_binom = binomtest(real_zero_count, real_total, p_null, alternative='two-sided')
    return real_zero_count, real_total, p_null, p_value_binom

def binomial_test_in_range(real, shuffled, lower=-10, upper=10, alternative='two-sided'):
    """
    Test whether the fraction of distances in [lower,upper] in `real` 
    differs from the null fraction in `shuffled`.
    
    Returns:
      real_count, real_n, p_null, p_value
    """
    # count real
    real_in = sum(1 for d in real if lower <= d <= upper)
    real_n  = len(real)
    # count null
    null_in = sum(1 for d in shuffled if lower <= d <= upper)
    null_n  = len(shuffled)
    p_null  = null_in / null_n if null_n>0 else 0.0

    bt = binomtest(real_in, real_n, p_null, alternative=alternative)
    return real_in, real_n, p_null, bt.pvalue

def compare_distances(real_distances, shuffled_distances_list):
    import numpy as np
    import scipy.stats as stats
    from statsmodels.stats.contingency_tables import Table2x2
    from scipy.stats import norm, ks_2samp
    from scipy.stats import mannwhitneyu, binomtest
    # Compute the mean distance for the real data
    #real_mean_distance = np.mean(real_distances)

    # Compute the mean distance for each shuffled iteration
    #shuffled_mean_distances = np.array([np.mean(distances) for distances in shuffled_distances_list])
    
    # Compute null distribution parameters (mean and std) from the shuffled data alone
    #mu_null = np.mean(shuffled_mean_distances)
    #sigma_null = np.std(shuffled_mean_distances)

    # Compute the Z-score relative to the null distribution
    #z_score_permutation = (real_mean_distance - mu_null) / sigma_null
    #p_value_permutation = 2 * norm.sf(abs(z_score_permutation))

    #extreme_count = sum(1 for d in shuffled_mean_distances if abs(d) >= abs(real_mean_distance))
    #p_value_permutation_extreme = (extreme_count + 1) / (len(shuffled_mean_distances) + 1)

    # Convert real_distances into 1D array
    real_distances_array = np.array(real_distances)
    # Convert shuffled distances into single 1D array
    shuffled_distances_array = np.concatenate(shuffled_distances_list)

    # Perform KS test
    ks_test_stat, ks_p_value = ks_2samp(real_distances_array, shuffled_distances_array)
    
    # Calculate additional distribution statistics
    real_median, real_iqr, real_std = get_distribution_stats(real_distances_array)
    shuffled_median, shuffled_iqr, shuffled_std = get_distribution_stats(shuffled_distances_array)
    
    # Mann-Whitney U test on a restricted region (e.g., -10 to 10; 35-55 spacing)
    #mw_stat, mw_p_value = mann_whitney_restricted(real_distances_array, shuffled_distances_array, lower_bound=-10, upper_bound=10)
    
    # MW over entire distribution (you probably should just take a slice of a distribution for mann-whitney since it is a rank order test)
    mw_stat, mw_p_value = mannwhitneyu(real_distances_array, shuffled_distances_array, alternative='two-sided')

    # Binomial test at 0 deviation (or consider a tolerance)
    real_zero_count, real_total, p_null_binom, binom_p_value = binomial_test_at_zero(real_distances_array, shuffled_distances_array, tolerance=0)

    # Binomial test at +/- 10 range 
    real_range_count, real_total, p_null_range_binom, binom_p_value_range = binomial_test_in_range(real_distances_array, shuffled_distances_array, lower=-10, upper=10, alternative='two-sided')


    # Chi squared test comparing null and real distributions at spacing=45
    # Count “ideal” (zero) vs “non‑ideal”
    a_zero = np.sum(real_distances_array == 0)
    b_zero = len(real_distances_array) - a_zero
    c_zero = np.sum(shuffled_distances_array == 0)
    d_zero = len(shuffled_distances_array) - c_zero

    # Count ideal over a range (35-55) vs "non-ideal"
    mask_real = (real_distances_array >= -10) & (real_distances_array <= 10)
    mask_null = (shuffled_distances_array >= -10) & (shuffled_distances_array <= 10)

    a_range = np.sum(mask_real)
    b_range = len(real_distances_array) - a_range
    c_range = np.sum(mask_null)
    d_range = len(shuffled_distances_array) - c_range


    # Construct contingency table
    contingency_table_zero = np.array([[a_zero, b_zero],
                                       [c_zero, d_zero]])

    # Perform chi-squared test without Yates' correction
    chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero = stats.chi2_contingency(contingency_table_zero, correction=False)
    
    # Calculate odds ratio and 95% confidence interval
    table2x2_zero = Table2x2(contingency_table_zero)
    odds_ratio_zero = table2x2_zero.oddsratio
    ci_low_zero, ci_upp_zero = table2x2_zero.oddsratio_confint()
    #p_value_oddsratio = table2x2.test_oddsratio_null().pvalue
    p_value_oddsratio_zero = table2x2_zero.oddsratio_pvalue()


    # Construct contingency table
    contingency_table_range = np.array([[a_range, b_range],
                                       [c_range, d_range]])

    # Perform chi-squared test without Yates' correction
    chi2_stat_range, p_value_chi_range, dof_range, expected_range = stats.chi2_contingency(contingency_table_range, correction=False)
    
    # Calculate odds ratio and 95% confidence interval
    table2x2_range = Table2x2(contingency_table_range)
    odds_ratio_range = table2x2_range.oddsratio
    ci_low_range, ci_upp_range = table2x2_range.oddsratio_confint()
    #p_value_oddsratio = table2x2.test_oddsratio_null().pvalue
    p_value_oddsratio_range = table2x2_range.oddsratio_pvalue()


    # Print the basic comparisons
    print(f'**********\n')
    print(f'OTHER TEST CATEGORIES')    
    print(f'KS test statistic: {ks_test_stat}')
    print(f'KS test p-value: {ks_p_value}\n')
    
    # Print additional distribution statistics
    print(f'Real median: {real_median}, IQR: {real_iqr}, Std: {real_std}')
    print(f'Shuffled median: {shuffled_median}, IQR: {shuffled_iqr}, Std: {shuffled_std}\n')
    
    # Print Mann-Whitney U test results
    print(f'Mann-Whitney U statistic (entire range): {mw_stat}')
    print(f'Mann-Whitney U p-value (entire range): {mw_p_value}\n')
    
    # Print Binomial test results at 0 deviation
    print(f'Real count at 0 deviation for biomial test: {real_zero_count} out of {real_total}')
    print(f'Estimated p(0 deviation) in shuffled distribution for binomial test: {p_null_binom}')
    print(f'Binomial two sided p-value at 0 deviation: {binom_p_value}\n')

    print(f"Real count in [-10,10] for binomial test: {real_range_count}/{real_total}")
    print(f"Null distribution count in [-10,10] for binomial test: {p_null_range_binom}")
    print(f"Binomial two‑sided p‑value for ±10 window: {binom_p_value_range}\n")


    # Print chi squared results (right at 45)
    print(f'Contingency table (deviation = 0 codons from ideal):')
    print(f'   Real ideal (zero)  = {a_zero}, real non-ideal (zero)  = {b_zero}')
    print(f'   Null ideal (zero)  = {c_zero}, null non-ideal (zero)  = {d_zero}')
    print(f'Chi squared stat (zero) = {chi2_stat_zero}; p-value (zero) = {p_value_chi_zero}')
    print(f'DOF (zero) = {dof_zero}; expected counts = {expected_zero}')
    print(f'Confidence interval 95% (lower, zero test) = {ci_low_zero}; Confidence interval 95% (upper, zero test) =  {ci_upp_zero}')
    print(f'Chi Squared OR (deviaion = 0) = {odds_ratio_zero:.3f}, OR p‑value (deviaion = 0) = {p_value_oddsratio_zero:.3g}\n')

    # Print chi squared results (35-55)
    print(f'Contingency table (deviation ±10 codons codons from ieal):')
    print(f'   Real ideal (35-55)  = {a_range}, real non-ideal (35-55)  = {b_range}')
    print(f'   Null ideal (35-55)  = {c_range}, null non-ideal (35-55)  = {d_range}')
    print(f'Chi squared stat (35-55) = {chi2_stat_range}; p-value (35-55) = {p_value_chi_range}')
    print(f'DOF (35-55) = {dof_range}; expected counts = {expected_range}')
    print(f'Confidence interval 95% (lower, 35-55 test) = {ci_low_range}; Confidence interval 95% (upper, 35-55 test) =  {ci_upp_range}')
    print(f'Chi Squared OR (deviaion -10 to 10) = {odds_ratio_range}, p‑value (deviaion -10 to 10) = {p_value_oddsratio_range}\n')

    print(f'**********\n')


    return (ks_test_stat, ks_p_value, real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std, mw_stat, mw_p_value, real_zero_count, real_total, p_null_binom, binom_p_value,
           real_range_count, p_null_range_binom, binom_p_value_range, 
           a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
           a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range)


def run_new_statistical_tests(tmd_data, filtered_list, slipsite_list, base_seed, output_file):

    strict_slip_sites = ['AAAAAAA', 'AAAAAAT', 'AAAAAAC', 'AAATTTA', 'AAATTTT', 'AAATTTC',
                         'TTTAAAA', 'TTTAAAT', 'TTTAAAC', 'TTTTTTA', 'TTTTTTT', 'TTTTTTC',
                         'CCCAAAA', 'CCCAAAT', 'CCCAAAC', 'CCCTTTA', 'CCCTTTT', 'CCCTTTC',
                         'GGGAAAA', 'GGGAAAT', 'GGGAAAC', 'GGGTTTA', 'GGGTTTT', 'GGGTTTC']

    # Get the subset of transcripts that are canonical only and alternative only
    filtered_list_canonical_only, filtered_list_alternative_only = filter_non_canonical_transcripts(filtered_list)

    # Fetch heptamer frequencies for the transcriptome
    heptamer_frequencies_all, set_unique_heptamers_all, shuffled_list_of_heptamers_all = calculate_heptamer_frequencies_prf(filtered_list, tmd_data)
    heptamer_frequencies_canonical, set_unique_heptamers_canonical, shuffled_list_of_heptamers_canonical = calculate_heptamer_frequencies_prf(filtered_list_canonical_only, tmd_data)
    heptamer_frequencies_alternative, set_unique_heptamers_alternative, shuffled_list_of_heptamers_alternative = calculate_heptamer_frequencies_prf(filtered_list_alternative_only, tmd_data)

    # Convert heptamer sets into lists for use in chi squared test
    all_heptamers_all = list(set_unique_heptamers_all)
    all_heptamers_canonical = list(set_unique_heptamers_canonical)
    all_heptamers_alternative = list(set_unique_heptamers_alternative)
    
    # Calculate p_null for each combination of heptamer frequencies and slip site lists
    p_null_strict_all = calculate_p_null_slippery(heptamer_frequencies_all, strict_slip_sites)
    p_null_strict_canonical = calculate_p_null_slippery(heptamer_frequencies_canonical, strict_slip_sites)
    p_null_strict_alternative = calculate_p_null_slippery(heptamer_frequencies_alternative, strict_slip_sites)

    p_null_loose_all = calculate_p_null_slippery(heptamer_frequencies_all, slipsite_list)
    p_null_loose_canonical = calculate_p_null_slippery(heptamer_frequencies_canonical, slipsite_list)
    p_null_loose_alternative = calculate_p_null_slippery(heptamer_frequencies_alternative, slipsite_list)

    # Log for probabilities of slip sites in categories
    print(f'Probability of strict slippery heptamers in all transcripts (p_null_strict_all): {p_null_strict_all}')
    print(f'Probability of strict slippery heptamers in canonical transcripts (p_null_strict_canonical): {p_null_strict_canonical}')
    print(f'Probability of strict slippery heptamers in alternative transcripts (p_null_strict_alternative): {p_null_strict_alternative}')
    print(f'Probability of loose slippery heptamers in all transcripts (p_null_loose_all): {p_null_loose_all}')
    print(f'Probability of loose slippery heptamers in canonical transcripts (p_null_loose_canonical): {p_null_loose_canonical}')
    print(f'Probability of loose slippery heptamers in alternative transcripts (p_null_loose_alternative): {p_null_loose_alternative}')


    # Initialize a list to store results
    results = []

    # Initialize dictionaries to store distance distributions for TMD shuffling and random heptamer sampling tests
    real_distances_tmd = {}
    shuffled_distributions_tmd = {}

    real_distances_heptamer = {}
    sampled_distributions_heptamer = {}


    # Initialize dictionaries to store ideal distances
    num_ideal_distances_real_tmd = {}
    num_ideal_distances_shuffled_tmd = {}


    # TMD Shuffling Test
    for category, data in [#('strict_all', (filtered_list, tmd_data, strict_slip_sites, base_seed)),
                           #('strict_canonical', (filtered_list_canonical_only, tmd_data, strict_slip_sites, base_seed)),
                           #('strict_alternative', (filtered_list_alternative_only, tmd_data, strict_slip_sites, base_seed)),
                           #('loose_all', (filtered_list, tmd_data, slipsite_list, base_seed)),
                           ('loose_canonical', (filtered_list_canonical_only, tmd_data, slipsite_list, base_seed)),
                           #('loose_alternative', (filtered_list_alternative_only, tmd_data, slipsite_list, base_seed))
                           ]:
        print(f"Performing TMD shuffling test on category: {category}")
        (z_score_permutation_tmd, 
        p_value_permutation_tmd, 
        real_mean_distance, mu_null, sigma_null, 
        extreme_count, p_value_permutation_extreme, 
        ks_test_stat, ks_p_value, 
        real_median, real_iqr, real_std, shuffled_median, shuffled_iqr, shuffled_std,
        mw_stat, mw_p_value, 
        real_zero_count, real_total, p_null_binom, binom_p_value,
        real_range_count, p_null_range_binom, binom_p_value_range,
        num_ideal_distances_real_tm, num_ideal_distances_shuffled, real_distances_tm, shuffled_distances_list,
        a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
        a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range) = compare_tmd_shuffle_to_observed(
            *data, num_iterations=16384)

        results.append({
            'test_type': 'TMD Shuffling',
            'category': category,
            'z_score_permutation': z_score_permutation_tmd,
            'p_z': p_value_permutation_tmd,
            'real_mean_distance' : real_mean_distance,
            'mu_null' : mu_null,
            'sigma_null' : sigma_null,
            'extreme_count' : extreme_count,
            'p_perm' : p_value_permutation_extreme,
            'ks_test_stat' : ks_test_stat,
            'ks_p_value' : ks_p_value,
            'real_median' : real_median,
            'real_iqr' : real_iqr,
            'real_std' : real_std,
            'shuffled_tmd_or_sampled_hept_median' : shuffled_median,
            'shuffled_tmd_or_sampled_hept_iqr' : shuffled_iqr,
            'shuffled_tmd_or_sampled_hept_std' : shuffled_std,
            'mann_whitney_stat' : mw_stat,
            'mann_whitney_p_value' : mw_p_value,
            'real_count_at_zero_deviation' : real_zero_count,
            'real_count_any_deviation_total' : real_total,
            'estimated_p_null_zero_in_shuffled_or_sampled' : p_null_binom,
            'binomial_test_at_zero_deviation_p_value' : binom_p_value,
            'real_count_35-55' : real_range_count,
            'estimated_p_null_35-55_in_shuffled_or_sampled' : p_null_range_binom,
            'binomial_test_at_35-55_deviation_p_value' : binom_p_value_range,
            'a_zero' : a_zero,
            'b_zero' : b_zero,
            'c_zero' : c_zero,
            'd_zero' : d_zero,
            'chi2_stat_zero' : chi2_stat_zero,
            'p_value_chi_zero' : p_value_chi_zero,
            'dof_zero' : dof_zero,
            'expected_counts_zero' : expected_zero,
            'ci_low_zero' : ci_low_zero,
            'ci_upper_zero' : ci_upp_zero,
            'odds_ratio_zero' : odds_ratio_zero,
            'p_value_oddsratio_zero' : p_value_oddsratio_zero,
            'a_range' : a_range,
            'b_range' : b_range,
            'c_range' : c_range,
            'd_range' : d_range,
            'chi2_stat_range' : chi2_stat_range,
            'p_value_chi_range' : p_value_chi_range,
            'dof_range' : dof_range,
            'expected_counts_range' : expected_range,
            'ci_low_range' : ci_low_range,
            'ci_upper_range' : ci_upp_range,
            'odds_ratio_range' : odds_ratio_range,
            'p_value_oddsratio_range' : p_value_oddsratio_range
        })

        num_ideal_distances_real_tmd[category] = num_ideal_distances_real_tm
        num_ideal_distances_shuffled_tmd[category] = num_ideal_distances_shuffled

        real_distances_tmd[category] = real_distances_tm
        shuffled_distributions_tmd[category] = shuffled_distances_list


    # Initialize dictionaries to store ideal distances
    num_ideal_distances_real_heptamer = {}
    num_ideal_distances_sampled_heptamer = {}


    # Heptamer Sampling Test
    for category, data, shuffled_heptamer_list in [
            #('strict_all', (filtered_list, tmd_data, strict_slip_sites, base_seed), shuffled_list_of_heptamers_all),
            #('strict_canonical', (filtered_list_canonical_only, tmd_data, strict_slip_sites, base_seed), shuffled_list_of_heptamers_canonical),
            #('strict_alternative', (filtered_list_alternative_only, tmd_data, strict_slip_sites, base_seed), shuffled_list_of_heptamers_alternative),
            #('loose_all', (filtered_list, tmd_data, slipsite_list, base_seed), shuffled_list_of_heptamers_all),
            ('loose_canonical', (filtered_list_canonical_only, tmd_data, slipsite_list, base_seed), shuffled_list_of_heptamers_canonical),
            #('loose_alternative', (filtered_list_alternative_only, tmd_data, slipsite_list, base_seed), shuffled_list_of_heptamers_alternative)
        ]:
        print(f"Performing heptamer sampling test on category: {category}")
        (z_score_permutation_hept, 
        p_value_permutation_hept, 
        real_mean_distance_hept, mu_null_hept, sigma_null_hept, 
        extreme_count_hept, p_value_permutation_extreme_hept, 
        ks_test_stat_hept, ks_p_value_hept, 
        real_median_hept, real_iqr_hept, real_std_hept, shuffled_median_hept, shuffled_iqr_hept, shuffled_std_hept,
        mw_stat_hept, mw_p_value_hept, 
        real_zero_count_hept, real_total_hept, p_null_binom_hept, binom_p_value_hept,
        real_range_count, p_null_range_binom, binom_p_value_range,
        num_ideal_distances_real_hept, num_ideal_distances_sampled, real_distances_hept, sampled_heptamer_distances_list,
        a_zero, b_zero, c_zero, d_zero, chi2_stat_zero, p_value_chi_zero, dof_zero, expected_zero, ci_low_zero, ci_upp_zero, odds_ratio_zero, p_value_oddsratio_zero,
        a_range, b_range, c_range, d_range, chi2_stat_range, p_value_chi_range, dof_range, expected_range, ci_low_range, ci_upp_range, odds_ratio_range, p_value_oddsratio_range) = compare_heptamer_sampling_to_observed(
            *data, shuffled_heptamer_list, num_iterations=16384)
        
        results.append({
            'test_type': 'Random Heptamer Sampling',
            'category': category,
            'z_score_permutation': z_score_permutation_hept,
            'p_z': p_value_permutation_hept,
            'real_mean_distance' : real_mean_distance_hept,
            'mu_null' : mu_null_hept,
            'sigma_null' : sigma_null_hept,
            'extreme_count' : extreme_count_hept,
            'p_perm' : p_value_permutation_extreme_hept,
            'ks_test_stat' : ks_test_stat_hept,
            'ks_p_value' : ks_p_value_hept,
            'real_median' : real_median_hept,
            'real_iqr' : real_iqr_hept,
            'real_std' : real_std_hept,
            'shuffled_tmd_or_sampled_hept_median' : shuffled_median_hept,
            'shuffled_tmd_or_sampled_hept_iqr' : shuffled_iqr_hept,
            'shuffled_tmd_or_sampled_hept_std' : shuffled_std_hept,
            'mann_whitney_stat' : mw_stat_hept,
            'mann_whitney_p_value' : mw_p_value_hept,
            'real_count_at_zero_deviation' : real_zero_count_hept,
            'real_count_any_deviation_total' : real_total_hept,
            'estimated_p_null_zero_in_shuffled_or_sampled' : p_null_binom_hept,
            'binomial_test_at_zero_deviation_p_value' : binom_p_value_hept,
            'real_count_35-55' : real_range_count,
            'estimated_p_null_35-55_in_shuffled_or_sampled' : p_null_range_binom,
            'binomial_test_at_35-55_deviation_p_value' : binom_p_value_range,
            'a_zero' : a_zero,
            'b_zero' : b_zero,
            'c_zero' : c_zero,
            'd_zero' : d_zero,
            'chi2_stat_zero' : chi2_stat_zero,
            'p_value_chi_zero' : p_value_chi_zero,
            'dof_zero' : dof_zero,
            'expected_counts_zero' : expected_zero,
            'ci_low_zero' : ci_low_zero,
            'ci_upper_zero' : ci_upp_zero,
            'odds_ratio_zero' : odds_ratio_zero,
            'p_value_oddsratio_zero' : p_value_oddsratio_zero,
            'a_range' : a_range,
            'b_range' : b_range,
            'c_range' : c_range,
            'd_range' : d_range,
            'chi2_stat_range' : chi2_stat_range,
            'p_value_chi_range' : p_value_chi_range,
            'dof_range' : dof_range,
            'expected_counts_range' : expected_range,
            'ci_low_range' : ci_low_range,
            'ci_upper_range' : ci_upp_range,
            'odds_ratio_range' : odds_ratio_range,
            'p_value_oddsratio_range' : p_value_oddsratio_range



        })
        num_ideal_distances_real_heptamer[category] = num_ideal_distances_real_hept
        num_ideal_distances_sampled_heptamer[category] = num_ideal_distances_sampled

        real_distances_heptamer[category] = real_distances_hept
        sampled_distributions_heptamer[category] = sampled_heptamer_distances_list

    # Chi-Squared Test
    chi_squared_results = []

    for category, data, slip_sites, all_heptamers in [
            #('strict_all', filtered_list, strict_slip_sites, all_heptamers_all),
            #('strict_canonical', filtered_list_canonical_only, strict_slip_sites, all_heptamers_canonical),
            #('strict_alternative', filtered_list_alternative_only, strict_slip_sites, all_heptamers_alternative),
            #('loose_all', filtered_list, slipsite_list, all_heptamers_all),
            ('loose_canonical', filtered_list_canonical_only, slipsite_list, all_heptamers_canonical),
            #('loose_alternative', filtered_list_alternative_only, slipsite_list, all_heptamers_alternative)
        ]:
        counts = analyze_slipsite_to_tmd_distances_contingency_table(
            data, tmd_data, slip_sites, all_heptamers)
        
        stats_results = perform_statistical_analysis(counts)
        
        chi_squared_results.append({
            'test_type': 'Chi-Squared',
            'category': category,
            'p_value': stats_results['p_value'],
            'chi2_stat': stats_results['chi2_stat'],
            'degrees_of_freedom': stats_results['dof'],
            'odds_ratio': stats_results['odds_ratio'],
            'ci_lower': stats_results['ci_lower'],
            'ci_upper': stats_results['ci_upper'],
            'A': counts['slippery_ideal'],       # Adding count A
            'B': counts['non_slippery_ideal'],   # Adding count B
            'C': counts['slippery_other'],       # Adding count C
            'D': counts['non_slippery_other']    # Adding count D
        })

    # Write to CSV
    fieldnames = ['test_type', 'category', 'p_value', 'observed_count', 'z_score', 'chi2_stat', 'degrees_of_freedom', 'odds_ratio', 'ci_lower', 'ci_upper', 'A', 'B', 'C', 'D']
    with open(output_file+'_chisquared_test_output_results.tsv', 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for result in chi_squared_results:
            # Fill missing fields with None
            for field in fieldnames:
                if field not in result:
                    result[field] = None
            writer.writerow(result)


    print(f"Chi squared test results written to {output_file+'_chisquared_test_output_results.tsv'}")
    # Additional logging or saving of intermediate data can be added here


    # CHANGE THIS IF YOU WANT TO WRITE OTHER CATEGORIES
    categories = [
                  #'strict_all',
                  #'strict_canonical',
                  #'strict_alternative', 
                  #'loose_all',
                  'loose_canonical',
                  #'loose_alternative'
                   ]

    # Determine the number of iterations
    num_iterations_tmd_shuffle = len(next(iter(num_ideal_distances_shuffled_tmd.values())))

    # Write to CSV
    with open(output_file+'_tmd_shuffling_ideal_distances.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration'] + categories)
        for i in range(num_iterations_tmd_shuffle):
            row = [i + 1]
            for category in categories:
                row.append(num_ideal_distances_shuffled_tmd[category][i])
            writer.writerow(row)

    # Determine the number of iterations
    num_iterations_heptamer_sampling = len(next(iter(num_ideal_distances_sampled_heptamer.values())))

    # Write to CSV
    with open(output_file+'_heptamer_sampling_ideal_distances.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration'] + categories)
        for i in range(num_iterations_heptamer_sampling):
            row = [i + 1]
            for category in categories:
                row.append(num_ideal_distances_sampled_heptamer[category][i])
            writer.writerow(row)

    with open(output_file+'_ideal_distances_real_log.txt', 'w') as logfile:
        logfile.write('Ideal distances in unshuffled data (TMD Shuffling Test):\n')
        for category, count in num_ideal_distances_real_tmd.items():
            logfile.write(f"{category}: {count}\n")
        logfile.write('\nIdeal distances in unshuffled data (Heptamer Sampling Test):\n')
        for category, count in num_ideal_distances_real_heptamer.items():
            logfile.write(f"{category}: {count}\n")


    fieldnames = ['test_type',
                      'category',
                      'z_score_permutation',
                      'p_z',
                      'real_mean_distance',
                      'mu_null',
                      'sigma_null',
                      'extreme_count',
                      'p_perm',
                      'ks_test_stat',
                      'ks_p_value',
                      'real_median',
                      'real_iqr',
                      'real_std',
                      'shuffled_tmd_or_sampled_hept_median',
                      'shuffled_tmd_or_sampled_hept_iqr',
                      'shuffled_tmd_or_sampled_hept_std',
                      'mann_whitney_stat',
                      'mann_whitney_p_value',
                      'real_count_at_zero_deviation',
                      'real_count_any_deviation_total',
                      'estimated_p_null_zero_in_shuffled_or_sampled',
                      'binomial_test_at_zero_deviation_p_value',
                      'real_count_35-55',
                      'estimated_p_null_35-55_in_shuffled_or_sampled',
                      'binomial_test_at_35-55_deviation_p_value',
                      'a_zero',
                      'b_zero',
                      'c_zero',
                      'd_zero',
                      'chi2_stat_zero',
                      'p_value_chi_zero',
                      'dof_zero',
                      'expected_counts_zero',
                      'ci_low_zero',
                      'ci_upper_zero',
                      'odds_ratio_zero',
                      'p_value_oddsratio_zero',
                      'a_range',
                      'b_range',
                      'c_range',
                      'd_range',
                      'chi2_stat_range',
                      'p_value_chi_range',
                      'dof_range',
                      'expected_counts_range',
                      'ci_low_range',
                      'ci_upper_range',
                      'odds_ratio_range',
                      'p_value_oddsratio_range']


    # build a DataFrame, force that exact column order, then write
    df = pd.DataFrame(results)
    df = df.reindex(columns=fieldnames)

    out_xlsx = output_file + '_all_tests_output_by_category.xlsx'
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        df.to_excel(writer,
                sheet_name='All_Statistical_Tests',
                index=False)
    print(f"Wrote all statistical tests to {out_xlsx}")


    with open(output_file+'_other_new_statistical_tests_categories.tsv', 'w', newline='') as tsvfile:
        fieldnames = ['test_type',
                      'category',
                      'z_score_permutation',
                      'p_z',
                      'real_mean_distance',
                      'mu_null',
                      'sigma_null',
                      'extreme_count',
                      'p_perm',
                      'ks_test_stat',
                      'ks_p_value',
                      'real_median',
                      'real_iqr',
                      'real_std',
                      'shuffled_tmd_or_sampled_hept_median',
                      'shuffled_tmd_or_sampled_hept_iqr',
                      'shuffled_tmd_or_sampled_hept_std',
                      'mann_whitney_stat',
                      'mann_whitney_p_value',
                      'real_count_at_zero_deviation',
                      'real_count_any_deviation_total',
                      'estimated_p_null_zero_in_shuffled_or_sampled',
                      'binomial_test_at_zero_deviation_p_value',
                      'real_count_35-55',
                      'estimated_p_null_35-55_in_shuffled_or_sampled',
                      'binomial_test_at_35-55_deviation_p_value',
                      'a_zero',
                      'b_zero',
                      'c_zero',
                      'd_zero',
                      'chi2_stat_zero',
                      'p_value_chi_zero',
                      'dof_zero',
                      'expected_counts_zero',
                      'ci_low_zero',
                      'ci_upper_zero',
                      'odds_ratio_zero',
                      'p_value_oddsratio_zero',
                      'a_range',
                      'b_range',
                      'c_range',
                      'd_range',
                      'chi2_stat_range',
                      'p_value_chi_range',
                      'dof_range',
                      'expected_counts_range',
                      'ci_low_range',
                      'ci_upper_range',
                      'odds_ratio_range',
                      'p_value_oddsratio_range']

        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for result in results:
            writer.writerow(result)

    write_distance_histograms(
        real_distances_tmd,
        shuffled_distributions_tmd,
        output_file + '_tmd_shuffle_histogram.xlsx',
        -10, 10
    )

    write_distance_histograms(
        real_distances_heptamer,
        sampled_distributions_heptamer,
        output_file + '_heptamer_sampling_histogram.xlsx',
        -10, 10
    )


    # Sample and write shuffled TMD - slippery heptamer distances
    with open(output_file+'_tmd_shuffling_distance_distributions.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Category', 'Sample_ID', 'Distance'])

        for category, all_iters in shuffled_distributions_tmd.items():
            #real_size = len(real_distances_tmd[category])
            for i, iteration in enumerate(all_iters):
                #sampled = random.sample(iteration, real_size) if len(iteration) >= real_size else iteration
                #for d in sampled:
                for d in iteration:
                    writer.writerow([category, i + 1, d])


    # Sample and write random heptamer-TMD distances
    with open(output_file+'_heptamer_sampling_distance_distributions.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Category', 'Sample_ID', 'Distance'])

        for category, all_iters in sampled_distributions_heptamer.items():
            #real_size = len(real_distances_heptamer[category])
            for i, iteration in enumerate(all_iters):
                #sampled = random.sample(iteration, real_size) if len(iteration) >= real_size else iteration
                #for d in sampled:
                for d in iteration:
                    writer.writerow([category, i + 1, d])

    with open(output_file+'_real_slippery_tmd_distances.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Category', 'Distance'])
        for category, distances in real_distances_tmd.items():
            for d in distances:
                writer.writerow([category, d])

    print(f"Results written to {output_file}")
    print("Ideal distances data saved to CSV and log files.")


    

    #return shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative

def write_distance_histograms(real_dict, null_dict, out_xlsx,
                              window_min=-10, window_max=10):
    """
    real_dict:  { category: [d1, d2, …] }
    null_dict:  { category: [[d11,d12,…], [d21,d22,…], …] }
    out_xlsx:   path to write the Excel workbook
    """
    BINS = list(range(window_min, window_max + 1))

    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        for cat, real_list in real_dict.items():
            real_arr = np.array(real_list, dtype=int)

            # 1) real counts
            cnts_real   = Counter(real_arr)
            real_counts = [cnts_real[b] for b in BINS]

            # 2) null per-iteration counts
            null_iters = null_dict[cat]  # e.g. list of 16 384 lists
            # build a DataFrame: rows=iters, cols=bins
            null_counts_iter = pd.DataFrame(
                ( Counter(it) for it in null_iters ),
                columns=BINS
            ).fillna(0).astype(int)

            # 3) null mean/std per bin
            null_mean = null_counts_iter.mean(axis=0)
            null_std  = null_counts_iter.std(axis=0, ddof=1)

            # 4) assemble summary table
            summary_df = pd.DataFrame({
                'Distance':   BINS,
                'Real_Count': real_counts,
                'Null_Mean':  null_mean.values,
                'Null_STD':   null_std.values,
            })

            # 5) write to its own sheet (31-char max)
            sheet = cat[:31]
            summary_df.to_excel(writer, sheet_name=sheet, index=False)


def perform_statistical_analysis(counts):
    """
    Performs chi-squared test and calculates odds ratio with confidence intervals from counts.

    Args:
        counts (dict): Dictionary with counts in the contingency table.

    Returns:
        dict: Dictionary with chi-squared statistic, p-value, degrees of freedom, odds ratio, 95% confidence interval.
    """
    import numpy as np
    import scipy.stats as stats
    from statsmodels.stats.contingency_tables import Table2x2

    A = counts['slippery_ideal']
    B = counts['non_slippery_ideal']
    C = counts['slippery_other']
    D = counts['non_slippery_other']

    # Continuity correction to avoid divide-by-zero error
    if A == 0 or B == 0 or C == 0 or D == 0:
        A += 0.5
        B += 0.5
        C += 0.5
        D += 0.5

    # Construct contingency table
    contingency_table = np.array([[A, B],
                                  [C, D]])

    # Perform chi-squared test without Yates' correction
    chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table, correction=False)
    
    # Calculate odds ratio and 95% confidence interval
    table2x2 = Table2x2(contingency_table)
    odds_ratio = table2x2.oddsratio
    ci_low, ci_upp = table2x2.oddsratio_confint()
    #p_value_oddsratio = table2x2.test_oddsratio_null().pvalue
    p_value_oddsratio = table2x2.oddsratio_pvalue()

    results = {
        'chi2_stat': chi2_stat,
        'p_value': p_value,
        'dof': dof,
        'odds_ratio': odds_ratio,
        'ci_lower': ci_low,
        'ci_upper': ci_upp,
        'p_value_oddsratio': p_value_oddsratio
    }

    return results



def check_duplicate_columns(excel_path):
    df = pd.read_excel(excel_path)
    column_names = df.columns
    duplicate_columns = column_names[column_names.duplicated()].unique()
    print(f"Column names: ", column_names)
    if len(duplicate_columns) > 0:
        print(f"Duplicate columns found: {duplicate_columns}")
    else:
        print("No duplicate columns found.")


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
    
    # Return all values in the same order as interactive_configuration
    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys


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
    
    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys


def main():
    """
    This version can be used to be just the stats - no Harrington motif printout?
    """

    global friction_data, scoring_sys
    base_seed = 1022
    random.seed(base_seed)
    np.random.seed(base_seed)

    if len(sys.argv) < 2:
        print("Usage: python script.py <config_file.ini>")
        config_mode_flag = int(input("Do you ant to revert to interactive mode? [1] yes or [2] no"))
        if config_mode_flag == 1:
            path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys = interactive_configuration()
        else:
            sys.exit(1)
    
    config_file = sys.argv[1]

    path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, path_ensembl_features, gap_near, gap_far, outname, scoring_sys = configparser_function(config_file)

    fasta_list = parse_fasta(path_fasta)

    stop_filtered_fasta_list = stop_codon_filter(fasta_list)
    slipsite_list = import_text_file(path_slipsites)
    tmd_data = parse_tmd_csv(path_tmd_csv)
    friction_data = parse_friction_csv(path_friction_csv)


    ensembl_features_dict = parse_ensembl_uniprot_canon(path_ensembl_features)

    print("Processing transcripts and matching to Uniprot sequences...")
    filtered_list, canonical_count, alt_spliced_count = fasta_transcript_filter(stop_filtered_fasta_list, ensembl_features_dict)

    ensembl_uniprot_dict = process_transcripts_with_uniprot(filtered_list, blastp_results_path, tmd_data, outname)
    print_ensembl_uniprot_dict(ensembl_uniprot_dict, outname)


    print("Running new stats tests including Chi Squared test")
 
    run_new_statistical_tests(tmd_data, filtered_list, slipsite_list, base_seed, outname+'_new_stats_tests')

    print("Finished running bank of new tests!")


if __name__ == '__main__':

    main()
