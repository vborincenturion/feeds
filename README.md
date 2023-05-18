# FEEDS, Food wastE biopEptiDe claSsifier: from genome to function

FEEDS is a tool that takes a bacteria genome or predicted proteins of yeast genomes to classify the secreted protease profile. It also digests protein substrate sequences to predict peptides and classify the biopeptide sequences with a novel machine learning method. This tool is useful for identifying potential bioactive compounds in food and discovering novel applications for waste management.

## Install

Create a conda environment with Python 3.9 and packages:

```conda create -n feeds python=3.9 pandas biopython scikit-learn```

Enter the environment:

```conda activate feeds```

Then, install the required tools/packages:

- Prodigal - ```conda install -c bioconda prodigal```
- Diamond - ```conda install -c bioconda diamond```
- RapidPeptideGeneration - ```pip3 install rpg```
- Tensorflow - ```pip3 install tensorflow-addons[tensorflow-cpu]```

*If you encounter difficulties while installing all the packages, an alternative approach would be to utilize "mamba" instead of the "conda" function.
To learn how to install "mamba", you can refer to the official documentation available at https://mamba.readthedocs.io/en/latest/installation.html.

Clone this repository:

```git clone https://github.com/vborincenturion/feeds.git```

Locate and change the "rpg_user.py" of your RapidPeptideGenerator tool for the one in the db folder of this github.

Finally, use the script "download_db.py" to download all the databases to run the feeds tool.
- Make sure to have the ".py" files usage for all users: 
   
   ```chmod 777 download_db.py```
   
   ```chmod 777 feeds.py```

## Usage

    feeds.py [-h | --help]
    
    feeds.py [-t <num_threads>] [-k <kingdom>] [-f_length <filter_length>] [-f_mol <filter_weight>] [-d <digest>]

Options:

    -h, --help          Show this help message and exit.
    
    -t <num_threads>    Number of threads to use in the diamond tool [default: 1].
    
    -k <kingdom>        Kingdom of genome/protein file. "bacteria" will run the prodigal command to predict ORFs, "yeast" will skip this step. Required=True.
    
    -f_length <filter_length>         Filter peptide sequence by the length (aa number). Required=No.
    
    -f_mol <filter_weight>         Filter peptide sequence by the molecular weight (KDal). Required=No.
    
    -d <digest>         Digestion mode of RapidPeptideGenerator tool. ["s" or "c"]. Required=True.
    
Extension permitted
- bacteria: ".fasta", ".fna", ".fa"
- yeast: ".fasta", ".faa"
- substrate: ".fasta", ".faa"

To use the tool, place the genome you wish to test into the genome folder of the corresponding kingdom. Then, add the substrate protein sequence FASTA file into the substrate folder.

Run the feeds.py script with the required options. For example, to run the tool with 4 threads, a bacteria genome file, with a molecular weight filter and in sequential digestion mode, use:

```python feeds.py -t 4 -k bacteria -f_mol 4 -d s``` 

This will generate a table with the protease genome profile, an output fasta file with the predicted peptides digested by each secreted genome enzyme with a molecular protein weight filter of <= 4 KDal and a table with the biopeptide classification.

- Multiple genomes and substrate sequence files can be placed in the folders.
- To use the tool for Yeast genomes, users need to provide the predicted ORFs of the genome, as Prodigal only predicts ORFs for bacteria.

## Output

- merged/diamond/: MEROPS ORFs annotation with the RapidPeptideGenerator enzyme ID (RPG_CS column)
- peptide_length_ranges.csv: peptide length range information generated after the digestion mode
- peptide_weight_ranges.csv: peptide molecular weight range information generated after the digestion mode
- results/no_filtered/: all peptides sequence files generated after the digestion mode without length or weight filter
- results/filtered/: all peptides sequence files generated after the digestion mode with length or weight filter
- results/prediction/: results of machine learning class predictions
