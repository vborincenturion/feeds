#!/usr/bin/env python

# Imports
import os
import subprocess
import re
import pandas as pd
import glob
import argparse
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def parse_args():
    parser = argparse.ArgumentParser(
        description='FEEDS, Food wastE biopEptiDe claSsifier: from genome to function.',
        epilog='Use the link https://github.com/vborincenturion/feeds for more information',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-k', '--kingdom', choices=['bacteria', 'yeast'], required=True,
                        help='Genome file kingdom (choose "bacteria" or "yeast")')
    parser.add_argument('-f_length', '--filter_length', type=int,
                    help='Filter peptide sequence with a minimum length')
    parser.add_argument('-f_mol', '--filter_mol', type=int,
                    help='Filter peptide sequence with a minimum molecular weight')
    parser.add_argument('-d', '--digest', choices=['s', 'c'], required=True,
                        help='"s" sequential mode and "c" concurrent mode of RapidPeptideGenerator tool')
    return parser.parse_args()


# Define the required directories
directories = ['diamond', 'orf_prediction', 'merged', 'merged/diamond', 'peptide', 'filtered', 'filtered/length', 'filtered/weight']

# Create the directories if they do not exist
for directory in directories:
    if not os.path.exists(directory):
        os.mkdir(directory)

args = parse_args()
threads = args.threads

kingdom = args.kingdom

filter_length = args.filter_length

filter_mol = args.filter_mol

digest = args.digest

#Prodigal
if kingdom == 'bacteria':
    input_dir = f"genome/{kingdom}"
    output_dir = "orf_prediction"
    file_extension = (".fasta", ".fna", ".fa")

    for file_name in os.listdir(input_dir):
        if file_name.endswith(file_extension):
            input_file = os.path.join(input_dir, file_name)
            output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}.faa")
            cmd = f"prodigal -i {input_file} -a {output_file}"
            os.system(cmd)

    #Diamond
    input_dir = "orf_prediction"
    output_dir = "diamond"
    database_path = "db/merops"
    file_extension = (".faa", ".fasta")
    num_threads = args.threads

    for file_name in os.listdir(input_dir):
        if file_name.endswith(file_extension):
            input_file = os.path.join(input_dir, file_name)
            output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}.csv")
            cmd = ["diamond", "blastp", "--more-sensitive", "-k", "1", "-f", "6", "qseqid", "sseqid", "pident", "--id", "90", "--query-cover", "85", "--subject-cover", "85", "-d", database_path, "-q", input_file, "-o", output_file, "-p", str(num_threads)]
            subprocess.run(cmd)
            # add column names to the output file
            subprocess.run(f"echo 'qseqid\tsseqid\tpident' | cat - {output_file} > temp && mv temp {output_file}", shell=True)

elif kingdom == 'yeast':
    input_dir = f"genome/{kingdom}"
    output_dir = "diamond"
    database_path = "db/merops"
    file_extension = ".faa"
    num_threads = args.threads

    for file_name in os.listdir(input_dir):
        if file_name.endswith(file_extension):
            input_file = os.path.join(input_dir, file_name)
            output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}.csv")
            cmd = ["diamond", "blastp", "--more-sensitive", "-k", "1", "-f", "6", "qseqid", "sseqid", "pident", "--id", "90", "--query-cover", "85", "--subject-cover", "85", "-d", database_path, "-q", input_file, "-o", output_file, "-p", str(num_threads)]
            subprocess.run(cmd)
            # add column names to the output file
            subprocess.run(f"echo 'qseqid\tsseqid\tpident' | cat - {output_file} > temp && mv temp {output_file}", shell=True)

else:
    print(f"Invalid kingdom: {kingdom}. Please choose 'bacteria' or 'yeast'")

# Set the directory containing the diamond output files
diamond_dir = "diamond"

# Set the filename pattern to match for diamond output files
diamond_pattern = "*.csv"

# Set the filename for the MEROPS data file
merops_file = "db/merops_database.tsv"

# Loop over the diamond output files and merge them with the MEROPS data
for diamond_file in glob.glob(diamond_dir + "/" + diamond_pattern):
    diamond_df = pd.read_csv(diamond_file, delimiter='\t')
    merops_df = pd.read_csv(merops_file, delimiter='\t')
    merged_df = pd.merge(diamond_df, merops_df, left_on=diamond_df.columns[1], right_on=merops_df.columns[1], how='left')
    merged_file = diamond_file.replace(".csv", "_merged.csv")
    merged_df.to_csv("merged/" + merged_file, index=False)

def extract_rpg_cs_values(file_path):
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path, delimiter=',')
    
    # Extract the unique and numeric characters from the RPG_CS column
    rpg_cs_values = set(df['RPG_CS'].dropna().unique()) - {"no-secrete", "no-info", "peptidase-inhibitor"}
    print("List of RPG_CS values:", rpg_cs_values)
    
    # Return the result as a set of characters
    return sorted(list(rpg_cs_values))
    
directory_path = 'merged/diamond/'

## Loop over all files in the directory
if digest == 'c':
    for file_name in os.listdir(directory_path):
        # Check if the file has the `_merged.csv` extension
        if file_name.endswith('_merged.csv'):
            # Construct the full file path
            file_path = os.path.join(directory_path, file_name)

            # Extract the RPG_CS values (excluding "no-secrete" and "no-info" and NA or nan values) for this file
            rpg_cs_values = extract_rpg_cs_values(file_path)

            # Loop over all fasta files in the `substrate` directory
            for substrate_file_name in os.listdir('substrate'):
                if substrate_file_name.endswith('.fasta') or substrate_file_name.endswith('.faa'):
                    # Generate the `rpg` command with the input and output file names and paths, the `-e` option with the RPG_CS values for this file, and the `-d` option with the actual database name
                    rpg_command = f"rpg -i substrate/{substrate_file_name} -o peptide/{file_name[:-4]}_{substrate_file_name[:-6]}_peptide.fasta -e {'    '.join(rpg_cs_values)} -d c"
                
                    # Run the `rpg` command using `os.system`
                    os.system(rpg_command)

elif digest == 's':
    for file_name in os.listdir(directory_path):
        # Check if the file has the `_merged.csv` extension
        if file_name.endswith('_merged.csv'):
            # Construct the full file path
            file_path = os.path.join(directory_path, file_name)

            # Extract the RPG_CS values (excluding "no-secrete" and "no-info" and NA or nan values) for this file
            rpg_cs_values = extract_rpg_cs_values(file_path)
            
            # Loop over all fasta files in the `substrate` directory
            for substrate_file_name in os.listdir('substrate'):
                if substrate_file_name.endswith('.fasta') or substrate_file_name.endswith('.faa'):
                    # Generate the `rpg` command with the input and output file names and paths, the `-e` option with the RPG_CS values for this file, and the `-d` option with the actual database name
                    rpg_command = f"rpg -i substrate/{substrate_file_name} -o peptide/{file_name[:-4]}_{substrate_file_name[:-6]}_peptide.fasta -e {'    '.join(rpg_cs_values)} -d s"
                
                    # Run the `rpg` command using `os.system`
                    os.system(rpg_command)   
                
# Path to peptide directory
peptide_dir = "peptide/"

# List all fasta files in the peptide directory
fasta_files = [file for file in os.listdir(peptide_dir) if file.endswith(".fasta")]

# Create a dictionary to store the count of length ranges for each fasta file
count_dict_l = {"Filename": [], "1-4 aa": [], "5-20 aa": [], "21-100 aa": [], "101-Max aa": [], "Total Sequences": [], "Max Length": []}

# Loop through all fasta files and count the length ranges
for fasta_file in fasta_files:
    with open(peptide_dir + fasta_file, "r") as f:
        sequence_lengths = []
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if sequence:
                    sequence_lengths.append(len(sequence))
                    sequence = ""
            else:
                sequence += line.strip()
        sequence_lengths.append(len(sequence))
        count_dict_l["Filename"].append(fasta_file)
        count_dict_l["1-4 aa"].append(sum(1 for length in sequence_lengths if length <= 4))
        count_dict_l["5-20 aa"].append(sum(1 for length in sequence_lengths if 5 <= length <= 20))
        count_dict_l["21-100 aa"].append(sum(1 for length in sequence_lengths if 21 <= length <= 100))
        count_dict_l["101-Max aa"].append(sum(1 for length in sequence_lengths if length > 100))
        count_dict_l["Total Sequences"].append(len(sequence_lengths))
        count_dict_l["Max Length"].append(max(sequence_lengths))

# Convert the dictionary to a pandas DataFrame
df_l = pd.DataFrame.from_dict(count_dict_l)

# Save the DataFrame to a table
df_l.to_csv("peptide_length_ranges.csv", index=False)

# Filter sequences by length
if filter_length:
    input_dir = "peptide"
    output_dir = "filtered/length"

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            with open(input_file, "r") as f, open(output_file, "w") as out:
                for record in SeqIO.parse(f, "fasta"):
                    if len(record.seq) < filter_length:
                        SeqIO.write(record, out, "fasta")
  
# Create a dictionary of amino acids and their corresponding molecular weights
amino_acid_weights = {
    'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15, 'E': 129.12,
    'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.20, 'F': 147.18, 'P': 97.12, 'S': 87.08, 'T': 101.11, 'W': 186.21,
    'Y': 163.18, 'V': 99.13
}

# Create a dictionary to store the count of length ranges for each fasta file
count_dict_w = {"Filename": [], "0-1 kDa": [], "1-3 kDa": [], "3-5 kDa": [], "5-10 kDa": [], "Total Sequences": [], "Max Weight": []}

# Loop through all fasta files and count the molecular weight ranges
for fasta_file in fasta_files:
    with open(peptide_dir + fasta_file, "r") as f:
        sequence_weights = []
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if sequence:
                    # Calculate the molecular weight of the sequence
                    pa = ProteinAnalysis(sequence)
                    sequence_weights.append(pa.molecular_weight())
                    sequence = ""
            else:
                sequence += line.strip()
        # Calculate the molecular weight of the last sequence
        pa = ProteinAnalysis(sequence)
        sequence_weights.append(pa.molecular_weight())
        count_dict_w["Filename"].append(fasta_file)
        count_dict_w["0-1 kDa"].append(sum(1 for weight in sequence_weights if weight <= 1000))
        count_dict_w["1-3 kDa"].append(sum(1 for weight in sequence_weights if 1000 < weight <= 3000))
        count_dict_w["3-5 kDa"].append(sum(1 for weight in sequence_weights if 3000 < weight <= 5000))
        count_dict_w["5-10 kDa"].append(sum(1 for weight in sequence_weights if 5000 < weight <= 10000))
        count_dict_w["Total Sequences"].append(len(sequence_weights))
        count_dict_w["Max Weight"].append(max(sequence_lengths))

# Convert the dictionary to a pandas DataFrame
df_w = pd.DataFrame.from_dict(count_dict_w)

# Save the DataFrame to a table
df_w.to_csv("peptide_weight_ranges.csv", index=False)

# Filter sequences by molecular weight
if filter_mol:
    input_dir = "peptide"
    output_dir = "filtered/weight"

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            with open(input_file, "r") as f, open(output_file, "w") as out:
                for record in SeqIO.parse(f, "fasta"):
                    mw = molecular_weight(record.seq, "protein") / 1000  # Calculate molecular weight in kDa
                    if mw < filter_mw:
                            SeqIO.write(record, out, "fasta")
