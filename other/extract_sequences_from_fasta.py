import pandas as pd
import numpy as np
import os
from Bio import SeqIO, SeqUtils
from sys import argv

s = argv[1]
sequences = pd.DataFrame()
folder = "path_to_folder/" + s + "/fasta_files/"

for file in os.listdir(folder):
    structure_name = file.split(".")[0]
    fasta = SeqIO.parse(open(folder + file),'fasta')
    for seq in fasta:
        sequences.loc[structure_name, seq.id] = seq.seq


sequences.to_csv("path_to_folder/" + s + "/" + s + "_all_sequences_from_fasta.csv")
