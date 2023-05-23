#!/usr/bin/env python

# Imports
import os
import subprocess
import re
import pandas as pd
import numpy as np
import glob
import argparse
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tensorflow.keras.models import load_model
import tensorflow_addons as tfa
from itertools import product
import pickle

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
    parser.add_argument('-f_mol', '--filter_mol', type=float,
                    help='Filter peptide sequence with a minimum molecular weight')
    parser.add_argument('-d', '--digest', choices=['s', 'c'], required=True,
                        help='"s" sequential mode and "c" concurrent mode of RapidPeptideGenerator tool')
    return parser.parse_args()


# Define the required directories
directories = ['diamond', 'orf_prediction', 'merged', 'merged/diamond', 'results', 'results/filtered', 'results/no_filtered', 'results/prediction']

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
            output_file = os.path.join(output_dir, f"{file_name[:-len(file_extension)]}faa")
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
                    rpg_command = f"rpg -i substrate/{substrate_file_name} -o results/no_filtered/{file_name[:-4]}_{substrate_file_name[:-6]}_peptide.fasta -e {'    '.join(rpg_cs_values)} -d c"
                
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
                    rpg_command = f"rpg -i substrate/{substrate_file_name} -o results/no_filtered/{file_name[:-4]}_{substrate_file_name[:-6]}_peptide.fasta -e {'    '.join(rpg_cs_values)} -d s"
                
                    # Run the `rpg` command using `os.system`
                    os.system(rpg_command)   
                
# Path to peptide directory
peptide_dir = "results/no_filtered/"

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
    input_dir = "results/no_filtered"
    output_dir = "results/filtered"

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            with open(input_file, "r") as f, open(output_file, "w") as out:
                for record in SeqIO.parse(f, "fasta"):
                    if len(record.seq) <= filter_length:
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
    input_dir = "results/no_filtered"
    output_dir = "results/filtered"

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            with open(input_file, "r") as f, open(output_file, "w") as out:
                for record in SeqIO.parse(f, "fasta"):
                    seq = str(record.seq)
                    mol_weight = sum(amino_acid_weights.get(aa, 0) for aa in seq)
                    kda_weight = mol_weight / 1000
                    if kda_weight <= filter_mol:
                        SeqIO.write(record, out, "fasta")

# Delete files with no sequences
output_dir = "results/filtered"

for filename in os.listdir(output_dir):
    file_path = os.path.join(output_dir, filename)
    if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
        os.remove(file_path)                        
                    
#amino acid position for sparse encoding
aa_pos = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19,"-":20}

#blosum62 substitution matrix
aa_blosum62 = {"A":np.asarray([4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,-100]),
          "C":np.asarray([0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,-100]),
          "D":np.asarray([-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,-100]),
          "E":np.asarray([-1,-4,-2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,-100]),
          "F":np.asarray([-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,-100]),
          "G":np.asarray([0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,-100]),
          "H":np.asarray([-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,-100]),
          "I":np.asarray([-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,-100]),
          "K":np.asarray([-1,-3,-1,-1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,-100]),
          "L":np.asarray([-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,-100]),
          "M":np.asarray([-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,-100]),
          "N":np.asarray([-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,-100]),
          "P":np.asarray([-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,-100]),
          "Q":np.asarray([-1,-3,-0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,-100]),
          "R":np.asarray([-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,-100]),
          "S":np.asarray([1,-1,0,0,-2,0,-1,-2,0,-2,-1,-1,-1,0,-1,4,1,-2,-3,-2,-100]),
          "T":np.asarray([0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,100]),
          "V":np.asarray([0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,-100]),
          "W":np.asarray([-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,-100]),
          "Y":np.asarray([-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7,-100]),
          "-":np.asarray([-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,100,-100,-100,0])}

#amino acid combination for 3-mer encoding
aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
combinations_list = ["".join(comb) for comb in product(aa, repeat=3)]

def threemers(peptide_list):  #threemers encoding function
    pept_vector = np.zeros((len(peptide_list),8000))
    counter_list = 0
    for peptide in peptide_list:
        counts = {aa_comb: 0 for aa_comb in combinations_list}
        for i in range(len(peptide)-2):
            window = peptide[i:i+3]
            if window in counts:
                counts[window] += 1
        pept_vector[counter_list]=list(counts.values())
        counter_list +=1
    return pept_vector

def denser_peptide(peptide_list): #dense encoding function
    pept_vector = np.zeros((len(peptide_list),100))
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i]=aa_pos[peptide[i]]
        counter_list +=1
    return pept_vector

def blosum_substitution(peptide_list): #blosum matrix encoding function
    pept_vector = np.zeros((len(peptide_list),100,21))
    counter = 0
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i]=aa_blosum62[peptide[i]]
            counter += 1
        for j in range(counter,100):
            pept_vector[counter_list,j]=aa_blosum62["-"]
        counter_list +=1
        counter = 0
    lista = np.zeros((len(peptide_list),2100))
    for i in range(len(pept_vector)):
        lista[i] = pept_vector[i].flatten()
    return lista

def sparse_peptide(peptide_list): #sparse peptide function
    pept_vector = np.zeros((len(peptide_list),100,21,1))
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i,aa_pos[peptide[i]]]=1
        for j in range(len(peptide),100):
            pept_vector[counter_list,j,aa_pos["-"]]=1
        counter_list +=1
    lista = np.zeros((len(peptide_list),2100))
    for i in range(len(pept_vector)):
        lista[i] = pept_vector[i].flatten()
    return lista

def vectorize_peptide(peptide_list): #sparse encoding for NN
    pept_vector = np.zeros((len(peptide_list),100,21,1))
    counter = 0
    counter_list = 0
    for peptide in peptide_list:
        if type(peptide) == float:
            continue
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i,aa_pos[peptide[i]]]=1
            counter += 1
        for j in range(counter,100):
            pept_vector[counter_list,j,aa_pos["-"]]=1
        counter_list +=1
        counter = 0
    counter = 0
    counter_list = 0
    if len(peptide_list) == 1:
        pept_vector = np.squeeze(pept_vector,axis=3)
    else:
        pept_vector = np.squeeze(pept_vector)
    return pept_vector

#load all the models(from sklearn and keras presaved models)
antidiabetic = pickle.load(open("db/Models/Antidiabetic.sav", 'rb'))
antihypertensive = pickle.load(open("db/Models/Antihypertensive.sav", 'rb'))
antimicrobial = load_model("db/Models/Antimicrobial.h5", compile = True)
antioxidant = pickle.load(open("db/Models/Antioxidant.sav", 'rb'))
cardiovascular = pickle.load(open("db/Models/Cardiovascular.sav", 'rb'))
celiac = pickle.load(open("db/Models/Celiac_disease.sav", 'rb'))
immunomodulatory = pickle.load(open("db/Models/Immunomodulatory.sav", 'rb'))
neuropeptides = pickle.load(open("db/Models/Neuropeptides.sav", 'rb'))
opioid = load_model("db/Models/Opioid.h5", compile = True)

motifs = {"antidiabetic":["DPNI","DPNIE","NIE","PNI","PNIE","VDPN","VDPNI","VDPNIE","PVDP","PVDPN","PVDPNI","PVDPNIE"],
          "antihypertensive":["LHLP","NLHL","PLW","TMP","TTMP","HLPL","KTT","KTTM","KTTMP","NLHLP","PLIY","PLIYP"],
          "Antimicrobial":["KLL","KFG","FQW","FQWQ","FQWQR","FQWQRN","KFGK","QWQRN","RKV","WQRN"],
          "antioxidant":["HDH","EGS","GWNI","KVLPVPQK","LLPH","AMRL","AMRLT","DAH","DEQ","EPD","EQA","GGGA","GPPGPPGPP","HPE","LPHH","MRL","MRLT","NGRF","NRPC","PPGPPGPP","SKVLPVPQK","WNIP"],
          "cardiovascular":["ALRR","EAVD","EAVDA","KALR","QEAV","RGDS","RQEA","ALRRQ","EAVDAL","KALRR","KKALR","QEAVD","QEAVDA","RQEAV","RRQE","RRQEA"],
          "celiac":["QPF","QQPQ","QPFP","QQPF","QQPFP","QPQQP","PQQPQ","QPFPQ","FPQPQ","PQPQQ"],
          "immunomodulatory":["KQD","KQDK","RKQD","RKQDK","FNQL","FYF","KPFKF","NKPF","NKPFK","NKPFKF","PFKF","PFNQ","PFNQL","RNK","RNKP","RNKPF","RNKPFK","RNKPFKF","TKFY","VTKF","VTKFY"],
          "neuropeptides":["NFLR","NFLRF","RNFLR","RNFLRF","LRFG","LRLR","LRLRF","RLRF","FLRFG","PFVEPI"],
          "opioid":["PFGF","YPFG","YPFGF","YGGF","GFLR","GGFLR","YGGFL","YGGFLR","GGFM","YGGFM"]}

modelli = [antidiabetic, antihypertensive, antioxidant, cardiovascular, celiac, immunomodulatory, neuropeptides]
modelli_NN = [antimicrobial, opioid]
modelli_nomi = ["antidiabetic", "antihypertensive", "antioxidant", "cardiovascular", "celiac", "immunomodulatory", "neuropeptides"]
modelli_nomi_NN = ["antimicrobial", "opioid"]

for folder in os.listdir('results/filtered'):
    for filename in os.listdir('results/filtered/'):
        if filename.endswith('.fasta'):
            #prepare the dataframe to be used for storing the info
            input_file = 'results/filtered/' + filename
            records = list(SeqIO.parse(input_file, "fasta"))
            data_dict = {"Header": [], "Sequence": []}
            for i, record in enumerate(records):
                data_dict["Header"].append(record.description)
                data_dict["Sequence"].append(str(record.seq))
            df = pd.DataFrame.from_dict(data_dict)
            df["Similarity search"] = np.nan
            df["Motif search"] = np.nan
            appaiati = []
            motivati = []
            #search on the unfiltered database if there is a 100% correspondance
            full_db = pd.read_csv("db/Peptide_all.csv")
            for i in range(len(df["Sequence"])):
                correspondance = []
                motif_found = []
                for j in range(len(full_db["Sequence"])):
                    if df.iloc[i][1] == full_db.iloc[j][1]:
                        correspondance.append(full_db.iloc[j][0])
                        continue
                for k in motifs:
                    if df.iloc[i][1] in motifs[k]:
                        motif_found = k
                if len(correspondance) == 0:
                    correspondance.append("")
                if len(motif_found) == 0:
                    motif_found.append("")
                motivati.append(motif_found)
                appaiati.append(correspondance)
            df["Similarity search"] = appaiati
            df["Motif search"] = motivati

            #encoding the sequences in different ways
            sparse_encoding = sparse_peptide(df["Sequence"])
            threemer_encoded = threemers(df["Sequence"])
            dense_encoding = denser_peptide(df["Sequence"])
            blosum_encoding = blosum_substitution(df["Sequence"])
            sparse_NN = vectorize_peptide(df["Sequence"])
            #predict using keras NN models
            for m in range(len(modelli_NN)):
                predizione = modelli_NN[m].predict(sparse_NN) 
                df[modelli_nomi_NN[m]] = predizione.tolist() 
            #predict with sklearn
            for m in range(len(modelli)):
                print(modelli_nomi[m])
                if modelli_nomi[m] in ["antioxidant","cardiovascular","celiac","immunomodulatory"]:
                    predizione = modelli[m].predict_proba(threemer_encoded)
                elif modelli_nomi[m] == "antidiabetic":
                    predizione = modelli[m].predict_proba(sparse_encoding)
                elif modelli_nomi[m] == "neuropeptides":
                    predizione = modelli[m].predict_proba(dense_encoding)
                elif modelli_nomi[m] == "antihypertensive":
                    predizione = modelli[m].predict_proba(blosum_encoding)
                df[modelli_nomi[m]] = predizione.tolist() 
            #save file with the same name as file in input+_predicted_peptides
            df.to_csv("results/prediction/"+filename+"_predicted_peptides.csv", index = False)
