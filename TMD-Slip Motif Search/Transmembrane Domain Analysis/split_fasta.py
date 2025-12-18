#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def split_fasta(input_file, num_parts):
    with open(input_file, 'r') as infile:
        sequences = infile.read().split('>')[1:]  # Read and split the file by '>'
    
    # Calculate the number of sequences per part
    total_sequences = len(sequences)
    sequences_per_part = total_sequences // num_parts
    remainder = total_sequences % num_parts
    
    part_sequences = []
    start = 0
    
    for i in range(num_parts):
        end = start + sequences_per_part + (1 if i < remainder else 0)
        part_sequences.append(sequences[start:end])
        start = end
    
    base_name = input_file.rsplit('.', 1)[0]
    
    for i, seqs in enumerate(part_sequences):
        output_file = f"{base_name}_part_{i+1:02}.fasta"
        with open(output_file, 'w') as outfile:
            for seq in seqs:
                outfile.write(f'>{seq}')
    
    print(f"Successfully split into {num_parts} parts.")

def main():
    parser = argparse.ArgumentParser(description="Split a FASTA file into multiple parts.")
    parser.add_argument("input_file", help="The input FASTA file")
    parser.add_argument("num_parts", type=int, help="The number of parts to split the file into")
    
    args = parser.parse_args()
    
    split_fasta(args.input_file, args.num_parts)

if __name__ == "__main__":
    main()
