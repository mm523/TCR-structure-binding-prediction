import pandas as pd
import numpy as np
from Bio.PDB import *
from Bio import SeqUtils

stcr = pd.read_csv("path_to_folder/pdb_original/stcr_summary_annotated_norepeats.csv", index_col = 0)
sequence_csv = pd.DataFrame()

for id in stcr.index.values:
    try:
        print(id)
        file = "path_to_folder/pdb_original/renumbered_pdb/" + id + "_renumbered.pdb"
        Parser = PDBParser()
        structure = Parser.get_structure(id, file)
        alpha = "A"
        beta = "B"
        peptide = "P"
        mhc1 = "M"

        chain_names = ["alpha", "beta", "peptide", "mhc"]

        f = open("path_to_folder/pdb_original/fasta_files/" + id + ".fasta", "w+")
        for i, chain in enumerate([alpha, beta, peptide, mhc1]):
            chain_seq = []
            chain_name = chain_names[i]
            f.write(">" + chain_name + "\n")
            for res in structure[0][chain]:
                if SeqUtils.seq1(res.resname) != "X":
                    chain_seq.append(SeqUtils.seq1(res.resname))
            f.write("".join(chain_seq) + "\n")
            sequence_csv.loc[id, chain_name] = "".join(chain_seq)
        f.close()
        print(sequence_csv)
    except Exception as e:
        print(e)
        continue

sequence_csv.to_csv("path_to_folder/pdb_original/pdb_original_all_sequences.csv")
