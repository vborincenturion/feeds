#!/usr/bin/env python

# Imports
import os
import subprocess
import re
import pandas as pd
import glob
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description='FEEDS, Food wastE biopEptiDe claSsifier: from genome to function.',
        epilog='Use the link https://github.com/vborincenturion/feeds for more information',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-k', '--kingdom', choices=['bacteria', 'yeast'], required=True,
                        help='Genome file kingdom (choose "bacteria" or "yeast")')
    parser.add_argument('-f', '--filter', choices=['yes', 'no'], required=True,
                        help='Filter peptide sequence with >20 aa')
    return parser.parse_args()


# Define the required directories
directories = ['diamond', 'orf_prediction', 'merged', 'merged/diamond', 'peptide', 'filtered']

# Create the directories if they do not exist
for directory in directories:
    if not os.path.exists(directory):
        os.mkdir(directory)

args = parse_args()
threads = args.threads

kingdom = args.kingdom

filter = args.filter

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
    rpg_cs_values = set(df['RPG_CS'].dropna().unique()) - {"no-secrete", "no-info"}
    
    # Return the result as a set of characters
    return sorted(list(rpg_cs_values))

directory_path = 'merged/diamond/'

## Loop over all files in the directory
for file_name in os.listdir(directory_path):
    # Check if the file has the `_merged.csv` extension
    if file_name.endswith('_merged.csv'):
        # Construct the full file path
        file_path = os.path.join(directory_path, file_name)

        # Extract the RPG_CS values (excluding "no-secrete" and "no-info" and NA or nan values) for this file
        rpg_cs_values = extract_rpg_cs_values(file_path)

        # Loop over all fasta files in the `substrate` directory
        for substrate_file_name in os.listdir('substrate'):
            if substrate_file_name.endswith('.fasta'):
                # Generate the `rpg` command with the input and output file names and paths, the `-e` option with the RPG_CS values for this file, and the `-d` option with the actual database name
                rpg_command = f"rpg -i substrate/{substrate_file_name} -o peptide/{file_name[:-4]}_{substrate_file_name[:-6]}_peptide.fasta -e {'    '.join(rpg_cs_values)} -d c"
                
                # Run the `rpg` command using `os.system`
                os.system(rpg_command)

# Path to peptide directory
peptide_dir = "peptide/"

# List all fasta files in the peptide directory
fasta_files = [file for file in os.listdir(peptide_dir) if file.endswith(".fasta")]

# Create a dictionary to store the count of length ranges for each fasta file
count_dict = {"Filename": [], "0-4": [], "5-20": [], "21-100": [], "101-Max": [], "Total Sequences": []}

# Loop through all fasta files and count the length ranges
for fasta_file in fasta_files:
    with open(peptide_dir + fasta_file, "r") as f:
        sequence_lengths = [len(line.strip()) for line in f.readlines()[1::2]]
        count_dict["Filename"].append(fasta_file)
        count_dict["0-4"].append(sum(1 for length in sequence_lengths if length <= 4))
        count_dict["5-20"].append(sum(1 for length in sequence_lengths if 5 <= length <= 20))
        count_dict["21-100"].append(sum(1 for length in sequence_lengths if 21 <= length <= 100))
        count_dict["101-max"].append(sum(1 for length in sequence_lengths if length > 101))
        count_dict["Total Sequences"].append(len(sequence_lengths))

# Convert the dictionary to a pandas DataFrame
df = pd.DataFrame.from_dict(count_dict)

# Save the DataFrame to a table
df.to_csv("peptide_length_ranges.csv", index=False)

# iterate through all the fasta files in peptide directory
if filter == 'yes':
    for filename in os.listdir('peptide'):
        if filename.endswith('.fasta'):
            with open(f'peptide/{filename}', 'r') as f_in, \
                    open(f'filtered/{filename}', 'w') as f_out:
                for record in SeqIO.parse(f_in, 'fasta'):
                    if len(record.seq) < 21:
                        SeqIO.write(record, f_out, 'fasta')
