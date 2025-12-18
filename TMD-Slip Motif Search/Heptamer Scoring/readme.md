Heptamer scoring (frameshift friction)

This folder contains the scripts used to assign thermodynamic frameshift friction scores (ΔG_FS) to nucleotide sequences. These scores are used in the manuscript to define slippery heptamers for the TMD-slip motif search.

Scoring is implemented in frameshift_routines_v10.py and uses Turner nearest-neighbor parameters via the companion module turner_energies.py. For each sequence, codon–anticodon interactions in the 0 frame and −1 frame are reconstructed at the ribosomal P site and A site, and the free-energy difference between these states is reported as the friction score.

Although the script supports multiple interactive analysis modes, only analysis option 3 is required to reproduce the manuscript results.

Running the heptamer scoring

Run the interactive script with:

python frameshift_routines_v10.py

When prompted, select analysis option [3] (upload text file containing list of heptamers and get friction scores for each sequence). Provide a path to a plain-text file containing one 7-nt heptamer per line, enter a job name, and select scoring system [2] (Turner nearest neighbor parameters).

Inputs

The only required input for analysis option [3] is a text file containing heptamer sequences (one per line). This may be an exhaustive list of all 16,384 possible 7-mers or a filtered subset.

Outputs

The script writes a CSV file named using the supplied job name. The table contains the input heptamer sequences and their corresponding friction scores. This file is used downstream as the path_friction_csv input to the motif-search pipeline.
