import numpy as np
import pandas as pd
from Bio.PDB import *
import os
from Bio import SeqUtils

atlas = pd.read_csv("./tcr_atlas/atlas_with_names.csv", index_col=0)

sequence_csv = pd.DataFrame(index = list(atlas.structure_id), columns = ["alpha", "beta", "peptide", "mhc"])

for i in atlas.index.values:
    if isinstance(atlas.loc[i, "template_PDB"], str):
        id = atlas.loc[i, "structure_id"]
        file = "./tcr_atlas/structures/designed_pdb/" + id + ".pdb"
        Parser = PDBParser()
        structure = Parser.get_structure(id, file)
        alpha = "D"
        beta = "E"
        peptide = "C"
        mhc1 = "A"

        chain_names = ["alpha", "beta", "peptide", "mhc"]

        f = open("./tcr_atlas/fasta_files/" + id + ".fasta", "w+")

        for i, chain in enumerate([alpha, beta, peptide, mhc1]):
            chain_seq = []
            chain_name = chain_names[i]
            f.write(">" + chain_name + "\n")
            for res in structure[0][chain]:
                if SeqUtils.seq1(res.resname) != "X":
                    chain_seq.append(SeqUtils.seq1(res.resname))
            f.write("".join(chain_seq) + "\n")
            sequence_csv.at[id, chain_name] = "".join(chain_seq)

        f.close()

sequence_csv.to_csv("./tcr_atlas/atlas_all_sequences.csv")
