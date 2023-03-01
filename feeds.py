#!/usr/bin/env python

# Imports
import os
import subprocess
import re
import pandas as pd
import glob
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='FEEDS: Food wastE biopEptiDe claSsifier: from genome to function')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-k', '--kingdom', choices=['bacteria', 'yeast'], required=True,
                        help='Genome file kingdom (choose "bacteria" or "yeast")')
    return parser.parse_args()


# Define the required directories
directories = ['diamond', 'orf_prediction', 'merged', 'merged/diamond', 'peptide']

# Create the directories if they do not exist
for directory in directories:
    if not os.path.exists(directory):
        os.mkdir(directory)

args = parse_args()
threads = args.threads
kingdom = args.kingdom

#Prodigal
if kingdom == 'bacteria':
    input_dir = f"genome/{kingdom}"
    output_dir = "orf_prediction"
    file_extension = ".fasta"

for file_name in os.listdir(input_dir):
    if file_name.endswith(file_extension):
        input_file = os.path.join(input_dir, file_name)
        output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}.faa")
        cmd = f"prodigal -i {input_file} -a {output_file}"
        os.system(cmd)

else:
    print(f"Skipping prodigal command for kingdom {kingdom}")
    
#Diamond
if kingdom == 'bacteria':
	input_dir = "orf_prediction"
	output_dir = "diamond"
	database_path = "db/merops"
	file_extension = ".faa"
	num_threads = args.threads
else:
	input_dir = "genome/{kingdom}"
	output_dir = "diamond"
	database_path = "db/merops"
	file_extension = ".faa"
	num_threads = args.threads

for file_name in os.listdir(input_dir):
    if file_name.endswith(".faa"):
        input_file = os.path.join(input_dir, file_name)
        output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}.csv")
        cmd = ["diamond", "blastp", "--more-sensitive", "-k", "1", "-f", "6", "qseqid", "sseqid", "pident", "--id", "90"
, "--query-cover", "85", "--subject-cover", "85", "-d", database_path, "-q", input_file, "-o", output_file]
        subprocess.run(cmd)
        # add column names to the output file
        subprocess.run(f"echo 'qseqid\tsseqid\tpident' | cat - {output_file} > temp && mv temp {output_file}", shell=True)

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

