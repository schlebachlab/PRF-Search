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
import multiprocessing as mp
import math
from scipy import stats
from scipy.stats import ks_2samp
from scipy.stats import binom_test
from scipy.stats import mannwhitneyu
from scipy.stats import poisson
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import groupby
import matplotlib.pyplot as plt
import frameshift_routines_v10 as frameshift
import time
import configparser
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
import difflib

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

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


def parse_motif_csv_notdone(path_friction_csv):
    friction_data = {}
    with open(path_friction_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            heptamer = row['sequence']
            friction_score = float(row['friction_score'])
            friction_data[heptamer] = friction_score
    return friction_data


def import_csv(path_csv):
    file_obj = open(path_csv, "r")
    data_raw = []
    for line in file_obj:
        line = line.rstrip()
        line = line.split(',')
        data_raw.append(line)
    return data_raw


def import_tsv(path_tsv):
    file_obj = open(path_tsv, "r")
    data_raw = []
    for line in file_obj:
        line = line.rstrip()
        line = line.split('\t')
        data_raw.append(line)
    return data_raw


def fetch_membrane_proteins_lists(uniprot_data_raw):
    key_word_membrane_list = []
    key_word_transmembrane_helix_list = []
    all_human_reviewed_list = []
    disease_association_list = []
    key_word_gpcr_list = []
    key_word_ion_channel_list = []

    uniprot_data_list = uniprot_data_raw[1:]
    for entry in uniprot_data_list:
        uniprot_accession_id = entry[0]
        uniprot_keywords_string = entry[12]
        uniprot_keywords_list = uniprot_keywords_string.split(';')
        
        all_human_reviewed_list.append(uniprot_accession_id)

        if 'Membrane' in uniprot_keywords_list:
            key_word_membrane_list.append(uniprot_accession_id)
        if 'Transmembrane helix' in uniprot_keywords_list:
            key_word_transmembrane_helix_list.append(uniprot_accession_id)
        if 'Disease variant' in uniprot_keywords_list:
            disease_association_list.append(uniprot_accession_id)
        if 'G-protein coupled receptor' in uniprot_keywords_list:
            key_word_gpcr_list.append(uniprot_accession_id)
        if 'Ion channel' in uniprot_keywords_list:
            key_word_ion_channel_list.append(uniprot_accession_id)


    return key_word_membrane_list, key_word_transmembrane_helix_list, disease_association_list, all_human_reviewed_list, key_word_gpcr_list, key_word_ion_channel_list


def parse_go_terms(uniprot_data_raw):
    # uniprot_headers = uniprot_data_raw[0]
    # 0 Entry
    # 1 Reviewed
    # 2 Entry Name
    # 3 Protein names
    # 4 Gene Names
    # 5 Organism
    # 6 Length
    # 7 Gene Ontology (biological process)
    # 8 Gene Ontology (cellular component)
    # 9 Involvement in disease
    # 10 Transmembrane
    # 11 Gene Ontology (molecular function)
    # 12 Keywords

    # Unique_go_terms
    go_term_description_dict_all = {}
    go_term_description_dict_biological_process = {}
    go_term_description_dict_cellular_component = {}
    go_term_description_dict_molecular_function = {}

    # Gene:GO term dictionaries
    gene_ontology_biological_process_dict = {}
    gene_ontology_cellular_component_dict = {}
    gene_ontology_molecular_function_dict = {}

    # Other output dictionaries
    uniprot_keywords_dict = {}
    uniprot_keywords_string_dict = {}
    disease_involvement_dict = {}
    protein_name_dict = {}

    # Loop to iterate through uniprot items
    uniprot_data_list = uniprot_data_raw[1:]
    for entry in uniprot_data_list:

        # Fetch basic terms from list
        uniprot_accession_id = entry[0]
        uniprot_entry_name = entry[2]
        descriptive_protein_name = entry[3]
        gene_names_string = entry[4]
        gene_names_list = gene_names_string.split(' ')

        # Define regex pattern to separate GO terms from their descriptions
        pattern = r'^(.*)\s+\[(GO:\d+)\]$'

        # Process biological process GO terms
        go_biological_process_string = entry[7]
        go_biological_process_raw_list = go_biological_process_string.split(';')
        go_biological_process_list = [x.strip() for x in go_biological_process_raw_list]
        go_ids_biological_process = []
        go_descriptions_biological_process = []
        for go_string in go_biological_process_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_biological_process.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_biological_process.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all GO terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictionary for biological process GO terms
                if go_term_id not in go_term_description_dict_biological_process:
                    go_term_description_dict_biological_process[go_term_id] = []
                    go_term_description_dict_biological_process[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_biological_process[go_term_id]:
                        go_term_description_dict_biological_process[go_term_id].append(go_description)



        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_biological_process_dict[uniprot_accession_id] = [go_biological_process_list, go_ids_biological_process, go_descriptions_biological_process]


        # Process cellular component GO terms
        go_cellular_component_string = entry[8]
        go_cellular_component_raw_list = go_cellular_component_string.split(';')
        go_cellular_component_list = [x.strip() for x in go_cellular_component_raw_list]
        go_ids_cellular_component = []
        go_descriptions_cellular_component = []
        for go_string in go_cellular_component_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_cellular_component.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_cellular_component.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictionary for cellular component GO terms
                if go_term_id not in go_term_description_dict_cellular_component:
                    go_term_description_dict_cellular_component[go_term_id] = []
                    go_term_description_dict_cellular_component[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_cellular_component[go_term_id]:
                        go_term_description_dict_cellular_component[go_term_id].append(go_description)
        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_cellular_component_dict[uniprot_accession_id] = [go_cellular_component_list, go_ids_cellular_component, go_descriptions_cellular_component]


        # Process molecular function GO terms
        go_molecular_function_string = entry[11]
        go_molecular_function_raw_list = go_molecular_function_string.split(';')
        go_molecular_function_list = [x.strip() for x in go_molecular_function_raw_list]
        go_ids_molecular_function = []
        go_descriptions_molecular_function = []
        for go_string in go_molecular_function_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_molecular_function.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_molecular_function.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictioanry for molecular function GO terms
                if go_term_id not in go_term_description_dict_molecular_function:
                    go_term_description_dict_molecular_function[go_term_id] = []
                    go_term_description_dict_molecular_function[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_molecular_function[go_term_id]:
                        go_term_description_dict_molecular_function[go_term_id].append(go_description)

        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_molecular_function_dict[uniprot_accession_id] = [go_molecular_function_list, go_ids_molecular_function, go_descriptions_molecular_function]


        # Process Uniprot Keywords
        uniprot_keywords_string = entry[12]
        uniprot_keywords_string_dict[uniprot_accession_id] = uniprot_keywords_string
        uniprot_keywords_list = uniprot_keywords_string.split(';')
        uniprot_keywords_dict[uniprot_accession_id] = uniprot_keywords_list

        # Process Uniprot annotations regarding disease involvement
        uniprot_disease_involvement_string = entry[9]
        uniprot_disease_involvement_list = uniprot_disease_involvement_string.split(';')

        if uniprot_disease_involvement_list[0] == '':
            num_associated_diseases = 'Number of associated diseases: 0'
            disease_involvement_output_entry = [num_associated_diseases, 'No disease involvement']
        else:
            num_associated_diseases = 'Number of associated diseases: ' + str(len(uniprot_disease_involvement_list))
            disease_involvement_output_entry = [num_associated_diseases]
            for entry in uniprot_disease_involvement_list:
                disease_involvement_output_entry.append(entry)

        disease_involvement_dict[uniprot_accession_id] = disease_involvement_output_entry

        # Process Uniprot names and protein descriptions
        protein_name_entry = [uniprot_entry_name, descriptive_protein_name, gene_names_string, gene_names_list]
        protein_name_dict[uniprot_accession_id] = protein_name_entry

    all_outputs = (go_term_description_dict_all,
                   go_term_description_dict_biological_process,
                   go_term_description_dict_cellular_component,
                   go_term_description_dict_molecular_function,
                   gene_ontology_biological_process_dict,
                   gene_ontology_cellular_component_dict,
                   gene_ontology_molecular_function_dict,
                   uniprot_keywords_dict,
                   disease_involvement_dict,
                   protein_name_dict,
                   uniprot_keywords_string_dict)
    
    return all_outputs


def fetch_uniprot_id_from_motif_subset(harrington_motif_subset):

    subset_uniprot_ids_list = []
    for entry in harrington_motif_subset:
        uniprot_id = entry[32]
        subset_uniprot_ids_list.append(uniprot_id)
    
    subset_uniprot_ids_set = set(subset_uniprot_ids_list)
    subset_unique_uniprot_ids_list = list(subset_uniprot_ids_set)

    return subset_unique_uniprot_ids_list

def count_unique_transcripts_harrington_motifs(harrington_motifs):
    unique_transcripts = []
    for entry in harrington_motifs:
        transcript_id = entry[0]
        if transcript_id not in unique_transcripts:
            unique_transcripts.append(transcript_id)
    num_unique_transcripts = len(unique_transcripts)
    return num_unique_transcripts

def count_alt_only_harrington_motifs(harrington_motifs):
    motif_gene_bins = {}
    for entry in harrington_motifs:
        gene_id = entry[1]
        if gene_id not in motif_gene_bins:
            motif_gene_bins[gene_id] = []
            motif_gene_bins[gene_id].append(entry)
        else:
            motif_gene_bins[gene_id].append(entry)
    
    unique_alt_motifs = []
    for gene_id, list_of_motifs in motif_gene_bins.items():

        unique_motifs_canonical_sets = []
        unique_motifs_alternative_sets = []
        alternative_motifs = []
        canonical_motifs = []
        for entry in list_of_motifs:
            canonical_status = entry[26]
            tmd_seq = entry[6]
            heptad_seq = entry[11]
            gap_size = entry[12]
            harrington_motif = [tmd_seq, heptad_seq, gap_size]
            harrington_motif_set = set(harrington_motif)
            if canonical_status == 'Canonical':
                canonical_motifs.append(entry)
                if harrington_motif_set not in unique_motifs_canonical_sets:
                    unique_motifs_canonical_sets.append(harrington_motif_set)
            else:
                alternative_motifs.append(entry)
                if harrington_motif_set not in unique_motifs_alternative_sets:
                    unique_motifs_alternative_sets.append(harrington_motif_set)

        for entry in alternative_motifs:
            tmd_seq = entry[6]
            heptad_seq = entry[11]
            gap_size = entry[12]
            harrington_motif = [tmd_seq, heptad_seq, gap_size]
            harrington_motif_set = set(harrington_motif)
            if harrington_motif_set not in unique_motifs_canonical_sets:
                unique_alt_motifs.append(entry)

    unique_alt_motifs_dict = {}
    for entry in unique_alt_motifs:
        key = frozenset(entry)
        value = 'motif_is_unique_to_alt_isoform'
        unique_alt_motifs_dict[key] = value

    return unique_alt_motifs, unique_alt_motifs_dict


def add_annotations_to_motifs(harrington_motifs, key_word_gpcr_list, key_word_ion_channel_list, disease_association_list, uniprot_keywords_string_dict, loose_alt_unique_motif_dict):
    # 'gpcr?', 'ion_channel?' 'uniprot_disease_involvement', 'unique_to_alt_isoform', 'uniprot_keywords'
    motifs_with_annotations = []

    for entry in harrington_motifs:
        motif_frozenset = frozenset(entry)
        uniprot_id = entry[32]
        if uniprot_id in key_word_gpcr_list:
            gpcr_status = 'gpcr'
        else:
            gpcr_status = '-'
        if uniprot_id in key_word_ion_channel_list:
            channel_status = 'ion_channel'
        else:
            channel_status = '-'
        if uniprot_id in disease_association_list:
            disease_status = 'disease_linked'
        else:
            disease_status = '-'
        if motif_frozenset in loose_alt_unique_motif_dict:
            unique_to_alt_isoform = 'motif_unique_to_alt_isoform'
        else:
            unique_to_alt_isoform = '-'
        if uniprot_id in uniprot_keywords_string_dict:
            uniprot_keywords_string = uniprot_keywords_string_dict[uniprot_id]
        else:
            uniprot_keywords_string = '-'
        
        new_elements = [gpcr_status, channel_status, disease_status, unique_to_alt_isoform, uniprot_keywords_string]
        new_entry = entry + new_elements
        motifs_with_annotations.append(new_entry)

    return motifs_with_annotations


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


def parse_ensembl_uniprot_dict(path_ensembl_uniprot_dict):
    ensembl_uniprot_dict = {}
    with open(path_ensembl_uniprot_dict, newline = '') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            transcript_id = row['ensembl_transcript_id']
            uniprot_accession_id = row['uniprot_accession_id']
            ensembl_uniprot_dict[transcript_id] = uniprot_accession_id
    return ensembl_uniprot_dict


def map_all_predicted_tmd_transcripts_to_uniprot_ids(path_tmds='topcons_ensembl_output_march_2025.csv', path_ensembl_uniprot_dict='motif_search_output_v73_motif_search_ensembl_uniprot_dict_full_blastp.csv'):
    tmd_data = parse_tmd_csv(path_tmds)
    uniprot_ensembl_dict = parse_ensembl_uniprot_dict(path_ensembl_uniprot_dict)

    uniprot_accession_ids_with_tmds = []

    for ensembl_transcript_id, uniprot_accession_id in uniprot_ensembl_dict.items():
        if ensembl_transcript_id in tmd_data:
            if uniprot_accession_id not in uniprot_accession_ids_with_tmds:
                uniprot_accession_ids_with_tmds.append(uniprot_accession_id)
    return uniprot_accession_ids_with_tmds


def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()


def get_go_counts(study_set, background_set, go_term_dict):
    """
    Get counts for GO term enrichment analysis for Fisher's Exact Test.
    """
    # GO term counts
    go_term_counts = {}
    
    # Iterate over background to count GO terms
    for uniprot_id in background_set:
        if uniprot_id in go_term_dict:
            for go_term in go_term_dict[uniprot_id][1]:  # Access GO term IDs from the dictionary
                if go_term not in go_term_counts:
                    go_term_counts[go_term] = {'a': 0, 'b': 0, 'c': 0, 'd': 0}
                if uniprot_id in study_set:
                    go_term_counts[go_term]['a'] += 1  # Study set with GO term
                else:
                    go_term_counts[go_term]['b'] += 1  # Background with GO term

    # Calculate c and d
    for go_term, counts in go_term_counts.items():
        counts['c'] = len(study_set) - counts['a']  # Study set without GO term
        counts['d'] = len(background_set) - len(study_set) - counts['b']  # Background without GO term

    return go_term_counts


def run_fisher_exact_test(go_term_counts, go_term_description_dict_all, go_category, gene_set):
    p_values = []
    go_terms = []
    log_odds_ratios = []
    go_descriptions = []
    go_categories = []
    gene_set_list = []
    for go_term, counts in go_term_counts.items():
        contingency_table = [[counts['a'], counts['b']],
                             [counts['c'], counts['d']]]
        odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
        if odds_ratio > 0:
            log_odds_ratio = math.log(odds_ratio)
        else:
            log_odds_ratio = None
        p_values.append(p_value)
        go_terms.append(go_term)
        log_odds_ratios.append(log_odds_ratio)
        # Here, we take the first description (if multiple) to keep it concise.
        desc_list = go_term_description_dict_all.get(go_term, ["NA"])
        go_description = desc_list[0] if desc_list else "NA"
        go_descriptions.append(go_description)
        go_categories.append(go_category)
        gene_set_list.append(gene_set)
    return go_terms, p_values, log_odds_ratios, go_descriptions, go_categories, gene_set_list

def run_fisher_exact_test_old(go_term_counts, go_term_description_dict_all, go_category, gene_set):
    """
    Perform Fisher's Exact Test for each GO term and collect p-values.
    """
    p_values = []
    go_terms = []
    log_odds_ratios = []
    go_descriptions = []
    go_categories = []
    gene_set_list = []
    
    for go_term, counts in go_term_counts.items():
        # Build the contingency table for Fisher's Exact Test
        contingency_table = [[counts['a'], counts['b']],
                             [counts['c'], counts['d']]]
        
        # Perform Fisher's Exact Test
        odds_ratio, p_value = stats.fisher_exact(contingency_table, alternative='greater')

        # Compute log odds ratio and safely handle cases where odds_ratio is zero or undefined
        if odds_ratio > 0:
            log_odds_ratio = math.log(odds_ratio)
        else:
            log_odds_ratio = None

        #log_odds_ratio = log_odds_ratio_prior # Is this what I want?

        p_values.append(p_value)
        go_terms.append(go_term)
        log_odds_ratios.append(log_odds_ratio)
        go_description = go_term_description_dict_all[go_term]
        go_descriptions.append(go_description)
        go_categories.append(go_category)
        gene_set_list.append(gene_set)
    
    return go_terms, p_values, log_odds_ratios, go_descriptions, go_categories, gene_set_list


def apply_fdr_correction(p_values, alpha=0.05):
    """
    Apply FDR correction using the Benjamini-Hochberg method.
    """
    _, p_adjusted, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    return p_adjusted



def gene_ontology_enrichment_function(motif_subset_uniprot_list, background_uniprot_list, gene_ontology_dict, go_term_description_dict_all, go_category, gene_set):
    """
    Run a Fisher's exact test and p-value FDR correction using a subset of Uniprot accession IDs associated
    with Harrington motifs, an appropriate set of background genes (membrane proteins, all reviewed, etc),
    and a dictionary associating Uniprot IDs with GO terms
    """

    # Get the GO term counts for Fisher's exact test
    go_term_counts = get_go_counts(motif_subset_uniprot_list, background_uniprot_list, gene_ontology_dict)

    # Run Fisher's exact test for category of GO terms and subset of Harrington motifs
    go_terms, p_values, log_odds_ratios, go_descriptions, go_categories, gene_set_list = run_fisher_exact_test(go_term_counts, go_term_description_dict_all, go_category, gene_set)

    # Apply FDR correction
    p_adjusted = apply_fdr_correction(p_values)

    # Pull out the counts in the contingency table for each term in the same order
    a_counts = [go_term_counts[t]['a'] for t in go_terms]
    b_counts = [go_term_counts[t]['b'] for t in go_terms]
    c_counts = [go_term_counts[t]['c'] for t in go_terms]
    d_counts = [go_term_counts[t]['d'] for t in go_terms]


    # Save the results as a list
    output_list = list(zip(go_terms, go_descriptions, log_odds_ratios, p_values, p_adjusted, go_categories, gene_set_list, a_counts, b_counts, c_counts, d_counts))

    return output_list

def keyword_enrichment_function(study_set, background_set, keyword, uniprot_keywords_dict, gene_set_label):
    a = sum(1 for uid in study_set if uid in uniprot_keywords_dict and keyword in uniprot_keywords_dict[uid])
    n = len(study_set)
    K = sum(1 for uid in background_set if uid in uniprot_keywords_dict and keyword in uniprot_keywords_dict[uid])
    N = len(background_set)
    table = [[a, n - a],
             [K - a, (N - n) - (K - a)]]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return (gene_set_label, odds_ratio, p_value)

def apply_fdr_correction_keyword(enrichment_results, alpha=0.05):
    p_values = [res[2] for res in enrichment_results]
    _, p_adjusted, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    corrected_results = []
    for i, res in enumerate(enrichment_results):
        corrected_results.append((res[0], res[1], res[2], p_adjusted[i]))
    return corrected_results


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

def print_motif_csv_with_annotations(output_array, outname):
    with open(outname + '.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Updated headers to match the structure of candidate_entry and include Uniprot information
        csv_writer.writerow([
            'transcript_id', 'gene_id', 'gene_name', 'tm_start', 'tm_end', 'tm_dg', 'tm_domain', 'tm_codons',
            'ss_coord_nuc', 'ss_coord_codon', 'chr_coordinate', 'slipsite_seq', 'gap_size', 'location_ideality',
            'gap_nuc', 'gap_%gc', 'gap_avg_codon_freq', 'fs_peptide', 'terminated', 'transcript_seq_nuc',
            'transcript_seq_res', 'friction_score', 'slipsite_canonical?', 'stem_loop_found', 'pseudoknot_found', 'secondary_structure',
            'canonical_status', 'tmd_type', 'tmd_orientation', 'ensembl_transcript_type', 'ensembl_transcript_id_aliases', 'ensembl_uniprot_id', 'blasp_uniprot_id', 'blasp_similarity_score',
            'gpcr?', 'ion_channel?' 'uniprot_disease_involvement', 'unique_to_alt_isoform', 'uniprot_keywords'
        ])
        
        # Write motif entries to the CSV
        for motif in output_array:
            csv_writer.writerow(motif)


def main_go_enrichment_function(path_motifs_csv, path_uniprot_tsv, blastp_results_path, path_tmd_csv, go_obo_path, outname):

    from goatools.obo_parser import GODag
    from goatools.go_enrichment import GOEnrichmentStudy

    godag = GODag(go_obo_path)

    # Load raw uniprot data
    uniprot_data_raw = import_tsv(path_uniprot_tsv)
    uniprot_data = uniprot_data_raw[1:]
    uniprot_go_dicts = parse_go_terms(uniprot_data_raw)

    # Unwrap GO dictionaries
    (go_term_description_dict_all, 
     go_term_description_dict_biological_process, 
     go_term_description_dict_cellular_component, 
     go_term_description_dict_molecular_function, 
     gene_ontology_biological_process_dict, 
     gene_ontology_cellular_component_dict, 
     gene_ontology_molecular_function_dict, 
     uniprot_keywords_dict, 
     disease_involvement_dict, 
     protein_name_dict,
     uniprot_keywords_string_dict) = uniprot_go_dicts

    # Build single gene->go-IDs mapping
    gene2go = {}
    for uni, (_, go_ids, _) in gene_ontology_biological_process_dict.items():
        gene2go.setdefault(uni, set()).update(go_ids)
    for uni, (_, go_ids, _) in gene_ontology_cellular_component_dict.items():
        gene2go.setdefault(uni, set()).update(go_ids)
    for uni, (_, go_ids, _) in gene_ontology_molecular_function_dict.items():
        gene2go.setdefault(uni, set()).update(go_ids)

    
    # Load Harrington motif data
    harrington_motifs_raw = import_csv(path_motifs_csv)
    harrington_motifs = harrington_motifs_raw[1:]

    [strict_all_motifs, 
     strict_canonical_motifs, 
     strict_alternative_motifs, 
     loose_all_motifs, 
     loose_canonical_motifs, 
     loose_alternative_motifs] = split_up_harrington_motifs(harrington_motifs, 20, False)
    
    # MOVE THIS TO OTHER ROUTINE TOO
    # Define loose and strict unqiue alternative motifs
    loose_alternative_unique_alt_motifs, loose_alt_unique_motif_dict = count_alt_only_harrington_motifs(loose_all_motifs)
    strict_alternative_unique_alt_motifs, strict_alt_unique_motif_dict = count_alt_only_harrington_motifs(strict_all_motifs)

    # Count transcripts from each class - move ot other script
    num_strict_all = count_unique_transcripts_harrington_motifs(strict_all_motifs)
    num_strict_canonical = count_unique_transcripts_harrington_motifs(strict_canonical_motifs)
    num_strict_alternative = count_unique_transcripts_harrington_motifs(strict_alternative_motifs)
    num_loose_all = count_unique_transcripts_harrington_motifs(loose_all_motifs)
    num_loose_canonical = count_unique_transcripts_harrington_motifs(loose_canonical_motifs)
    num_loose_alternative = count_unique_transcripts_harrington_motifs(loose_alternative_motifs)
    num_uniquely_alt_loose = count_unique_transcripts_harrington_motifs(loose_alternative_unique_alt_motifs)
    num_uniquely_alt_strict = count_unique_transcripts_harrington_motifs(strict_alternative_unique_alt_motifs)

    # Get lists of Uniprot accession IDs for background
    (key_word_membrane_list, 
     key_word_transmembrane_helix_list_not_used, # Change this variable name so it is not passed to gene_ontology_enrichment_function
     disease_association_list, 
     all_human_reviewed_list,
     key_word_gpcr_list,
     key_word_ion_channel_list) = fetch_membrane_proteins_lists(uniprot_data_raw)

    gpcr_disease_logfile = open(f'{outname}_motifs_spacing_only_disease_and_gpcr_log.txt', 'w')

    # Get Uniprot IDs of Harrington motif subsets
    strict_all_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_all_motifs)
    strict_canonical_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_canonical_motifs)
    strict_alternative_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_alternative_motifs)
    loose_all_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_all_motifs)
    loose_canonical_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_canonical_motifs)
    loose_alternative_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_alternative_motifs)

    loose_alternative_unique_alt_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_alternative_unique_alt_motifs)
    strict_alternative_unique_alt_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_alternative_unique_alt_motifs)

    # Find number of GPCR and disease associated harrington motifs for sets
    strict_all_disease_association_list = list(set(disease_association_list) & set(strict_all_motifs_uniprot_ids))
    strict_canonical_disease_association_list = list(set(disease_association_list) & set(strict_canonical_motifs_uniprot_ids))
    strict_alternative_disease_association_list = list(set(disease_association_list) & set(strict_alternative_motifs_uniprot_ids))
    loose_all_disease_association_list = list(set(disease_association_list) & set(loose_all_motifs_uniprot_ids))
    loose_canonical_disease_association_list = list(set(disease_association_list) & set(loose_canonical_motifs_uniprot_ids))
    loose_alternative_disease_association_list = list(set(disease_association_list) & set(loose_alternative_motifs_uniprot_ids))

    print(f'Number of genes within Uniprot containing the keyword \'Disease Variant\': {len(disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of unique genes within Uniprot containing the keyword \'Disease Variant\': {len(set(disease_association_list))}\n', file=gpcr_disease_logfile)

    strict_all_gpcr_list = list(set(key_word_gpcr_list) & set(strict_all_motifs_uniprot_ids))
    strict_canonical_gpcr_list = list(set(key_word_gpcr_list) & set(strict_canonical_motifs_uniprot_ids))
    strict_alternative_gpcr_list = list(set(key_word_gpcr_list) & set(strict_alternative_motifs_uniprot_ids))
    loose_all_gpcr_list = list(set(key_word_gpcr_list) & set(loose_all_motifs_uniprot_ids))
    loose_canonical_gpcr_list = list(set(key_word_gpcr_list) & set(loose_canonical_motifs_uniprot_ids))
    loose_alternative_gpcr_list = list(set(key_word_gpcr_list) & set(loose_alternative_motifs_uniprot_ids))

    print(f'Number of genes within Uniprot containing the keyword \'G-protein coupled receptor\': {len(key_word_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of unique genes within Uniprot containing the keyword \'G-protein coupled receptor\': {len(set(key_word_gpcr_list))}\n', file=gpcr_disease_logfile)

    strict_all_ionchannel_list = list(set(key_word_ion_channel_list) & set(strict_all_motifs_uniprot_ids))
    strict_canonical_ionchannel_list = list(set(key_word_ion_channel_list) & set(strict_canonical_motifs_uniprot_ids))
    strict_alternative_ionchannel_list = list(set(key_word_ion_channel_list) & set(strict_alternative_motifs_uniprot_ids))
    loose_all_ionchannel_list = list(set(key_word_ion_channel_list) & set(loose_all_motifs_uniprot_ids))
    loose_canonical_ionchannel_list = list(set(key_word_ion_channel_list) & set(loose_canonical_motifs_uniprot_ids))
    loose_alternative_ionchannel_list = list(set(key_word_ion_channel_list) & set(loose_alternative_motifs_uniprot_ids))

    print(f'Number of genes within Uniprot containing the keyword \'Ion Channel\': {len(key_word_ion_channel_list)}', file=gpcr_disease_logfile)
    print(f'Number of unique genes within Uniprot containing the keyword \'Ion Channel\': {len(set(key_word_ion_channel_list))}\n', file=gpcr_disease_logfile)


    # Print out size of lists
    print(f'Number of disease associated genes within strict_all motifs: {len(strict_all_disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of disease associated genes within strict_canonical motifs: {len(strict_canonical_disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of disease associated genes within strict_alternative motifs: {len(strict_alternative_disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of disease associated genes within loose_all motifs: {len(loose_all_disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of disease associated genes within loose_canonical motifs: {len(loose_canonical_disease_association_list)}', file=gpcr_disease_logfile)
    print(f'Number of disease associated genes within loose_alternative motifs: {len(loose_alternative_disease_association_list)}\n', file=gpcr_disease_logfile)

    print(f'Number of GPCRs within strict_all motifs: {len(strict_all_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of GPCRs within strict_canonical motifs: {len(strict_canonical_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of GPCRs within strict_alternative motifs: {len(strict_alternative_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of GPCRs within loose_all motifs: {len(loose_all_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of GPCRs within loose_canonical motifs: {len(loose_canonical_gpcr_list)}', file=gpcr_disease_logfile)
    print(f'Number of GPCRs within loose_alternative motifs: {len(loose_alternative_gpcr_list)}\n', file=gpcr_disease_logfile) 


    print(f'Number of ion channels within strict_all motifs: {len(strict_all_ionchannel_list)}', file=gpcr_disease_logfile)
    print(f'Number of ion channels within strict_canonical motifs: {len(strict_canonical_ionchannel_list)}', file=gpcr_disease_logfile)
    print(f'Number of ion channels within strict_alternative motifs: {len(strict_alternative_ionchannel_list)}', file=gpcr_disease_logfile)
    print(f'Number of ion channels within loose_all motifs: {len(loose_all_ionchannel_list)}', file=gpcr_disease_logfile)
    print(f'Number of ion channels within loose_canonical motifs: {len(loose_canonical_ionchannel_list)}', file=gpcr_disease_logfile)
    print(f'Number of ion channels within loose_alternative motifs: {len(loose_alternative_ionchannel_list)}\n', file=gpcr_disease_logfile)


    print(f'Number of motifs unique to alternative isoforms (loose slip sites): {len(loose_alternative_unique_alt_motifs)}', file=gpcr_disease_logfile)
    print(f'Number of motifs unique to alternative isoforms (strict slip sites):  {len(strict_alternative_unique_alt_motifs)}\n', file=gpcr_disease_logfile)

    print(f'Number of transcripts containing strict_all motifs: {num_strict_all}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing strict_canonical motifs: {num_strict_canonical}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing strict_alternative motifs: {num_strict_alternative}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing loose_all motifs: {num_loose_all}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing loose_canonical motifs: {num_loose_canonical}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing loose_alternative motifs: {num_loose_alternative}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing loose_uniquely_alternative motifs: {num_uniquely_alt_loose}', file=gpcr_disease_logfile)
    print(f'Number of transcripts containing strict_uniquely_alternative motifs: {num_uniquely_alt_strict}', file=gpcr_disease_logfile)
    
    gpcr_disease_logfile.close()

    motifs_with_annotations = add_annotations_to_motifs(harrington_motifs, key_word_gpcr_list, key_word_ion_channel_list, disease_association_list, uniprot_keywords_string_dict, loose_alt_unique_motif_dict)

    # Printing table of motifs unique to alternative isoforms
    print_motif_csv(loose_alternative_unique_alt_motifs, f'{outname}_motifs_unique_to_alternative_isoforms')
    print_motif_csv_with_annotations(motifs_with_annotations, f'{outname}_motif_table_with_new_annotations')

    # Re-assign background list of Uniprot sequences
    uniprot_accession_ids_with_tmds = map_all_predicted_tmd_transcripts_to_uniprot_ids(path_tmd_csv, blastp_results_path_csv)
    key_word_transmembrane_helix_list = uniprot_accession_ids_with_tmds

    # Initialize GOEA with background
    pop_ids = key_word_transmembrane_helix_list 
    goea = GOEnrichmentStudy(
        pop_ids,      # population / background
        gene2go,      # the full UniProt→GO mapping
        godag,        # your parsed ontology
        propagate_counts=True,      # count children as well
        alpha=0.05,                # significance cut-off for reporting
        methods=["fdr_bh"]         # multiple testing correction
    )   


    # Run Fisher's exact test on all six motif subsets 

    # Initialize a dictionary to hold DataFrames
    dataframes = {}

    # Define categories and GO categories
    categories = [
        ('Strict_All', strict_all_motifs_uniprot_ids),
        ('Strict_Canonical', strict_canonical_motifs_uniprot_ids),
        ('Strict_Alternative', strict_alternative_motifs_uniprot_ids),
        ('Strict_Uniquely_Alternative', strict_alternative_unique_alt_motifs_uniprot_ids),
        ('Loose_All', loose_all_motifs_uniprot_ids),
        ('Loose_Canonical', loose_canonical_motifs_uniprot_ids),
        ('Loose_Alternative', loose_alternative_motifs_uniprot_ids),
        ('Loose_Uniquely_Alternative', loose_alternative_unique_alt_motifs_uniprot_ids)
    ]

    go_categories = [
        ('Biological_Process', gene_ontology_biological_process_dict),
        ('Cellular_Component', gene_ontology_cellular_component_dict),
        ('Molecular_Function', gene_ontology_molecular_function_dict)
    ]

    '''
    # Loop over categories and GO categories to create DataFrames
    for gene_set_name, motif_uniprot_ids in categories:
        for go_category_name, go_dict in go_categories:
            output_list = gene_ontology_enrichment_function(
                motif_uniprot_ids, key_word_transmembrane_helix_list, go_dict, go_term_description_dict_all, go_category_name, gene_set_name
            )
            # Convert output_list to DataFrame
            df = pd.DataFrame(output_list, columns=['go_term', 'go_description', 'log_odds_ratio', 'p_value', 'p_adjusted', 'go_category', 'gene_set', 'go_term_present_counts_study_a', 'go_term_present_counts_background_only_b', 'go_term_absent_counts_study_c', 'go_term_absent_counts_background_only_d'])
            # Create a key for the dataframes dictionary
            key = f"{gene_set_name}_{go_category_name}"
            dataframes[key] = df
    '''

    all_results = {}
    for gene_set_name, study_ids in categories:
        res = goea.run_study(study_ids)
        for NS, label in [("BP","Biological_Process"),
                          ("CC","Cellular_Component"),
                          ("MF","Molecular_Function")]:
            rows = []
            for r in res:
                if r.NS != NS:
                    continue
                a = r.study_count
                b = r.pop_count - r.study_count
                c = r.study_n - r.study_count
                d = (r.pop_n - r.pop_count) - c


                contingency_table = [[a, b],
                                    [c, d]]
                oddsratio, p_value = fisher_exact(contingency_table, alternative='greater')
                if oddsratio > 0:
                    log_oddsratio = math.log(oddsratio)
                else:
                    log_oddsratio = None

                # recompute oddds ratio
                '''
                if b > 0 and c > 0:
                    oddsratio = (a * d) / (b * c)
                    log_oddsratio = math.log(oddsratio)
                else:
                    oddsratio = None
                    log_oddsratio = None
                '''
                
                rows.append({
                    "go_term" :                 r.GO,
                    "go_description" :          r.name,
                    #"definition" :              godag[r.GO].defn.text,
                    #"synonyms" :                ";".join(godag[r.GO].synonyms),
                    "go_category":              label,
                    "gene_set" :                gene_set_name,
                    "odds_ratio"    :           oddsratio,
                    "log_odds_ratio" :          log_oddsratio,
                    "fold_enrichment" :         (a / r.study_n) / (r.pop_count / r.pop_n) if r.pop_count and r.pop_n else None,
                    "ratio_in_study" :          a / r.study_n,
                    "ratio_in_pop" :            r.pop_count / r.pop_n,
                    "p_value" :                 r.p_uncorrected,
                    "p_adjusted" :              r.p_fdr_bh,
                    "study_count_A" :           a,
                    "pop_count_minus_study_count_B" : b,
                    "study_n_minus_study_count_C" :        c,
                    "pop_n_minus_pop_count_minus_c_D" :   d,
                    "background_count" :        r.pop_count,
                    "background_n"  :           r.pop_n,
                    "study_n"   :               r.study_n,
                    "GO_term_depth" :           godag[r.GO].depth,
                    "study_genes" :             ";".join(sorted(r.study_items))
                })
            all_results[f"{gene_set_name}_{label}"] = pd.DataFrame(rows)

    # create dataframes
    for sheet, df in all_results.items():
        dataframes[sheet] = df


    # Now perform keyword enrichment for the Uniprot keyword "Disease variant"
    keyword_results = []
    for gene_set_name, motif_uniprot_ids in categories:
        res = keyword_enrichment_function(motif_uniprot_ids, key_word_transmembrane_helix_list, 
                                           "Disease variant", uniprot_keywords_dict, gene_set_name)
        keyword_results.append(res)
    
    keyword_results_corrected = apply_fdr_correction_keyword(keyword_results)
    df_keyword = pd.DataFrame(keyword_results_corrected, columns=['gene_set', 'odds_ratio', 'p_value', 'p_adjusted'])
    dataframes['Disease_variant_keyword_enrichment'] = df_keyword



    # Write all DataFrames to an Excel file with separate sheets
    with pd.ExcelWriter(f'{outname}_gene_ontology_results.xlsx') as writer:
        for sheet_name, df in dataframes.items():
            # Ensure sheet names are less than 31 characters and unique
            sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

def configparser_function(config_file="config.ini"):
    """
    Loads configuration from an ini file using configparser.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Paths
    path_motifs_csv = config.get('Paths', 'path_motifs_csv')
    path_uniprot_tsv = config.get('Paths', 'path_uniprot_tsv')
    blastp_results_path_csv = config.get('Paths', 'blastp_results_path_csv')
    path_tmd_csv = config.get('Paths', 'path_tmd_csv')
    go_obo_path = config.get('Paths', 'gene_ontology_obo_path')
    
    # Parameters
    outname = config.get('Parameters', 'outname')
    
    # Return all values in the same order as interactive_configuration
    return path_motifs_csv, path_uniprot_tsv, blastp_results_path, path_tmd_csv, go_obo_path, outname 

def main():
    
    if len(sys.argv) < 2:
        print("Usage: python this_script.py <config_file.ini>")
        else:
            sys.exit(1)

    config_file = sys.argv[1]

    path_motifs_csv, path_uniprot_tsv, blastp_results_path, path_tmd_csv, go_obo_path, outname = configparser_function(config_file)
    
    main_go_enrichment_function(path_motifs_csv, path_uniprot_tsv, blastp_results_path, path_tmd_csv, go_obo_path, outname)    

if __name__ == '__main__':

    main()

