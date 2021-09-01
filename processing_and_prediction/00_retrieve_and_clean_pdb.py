from sys import argv
import pandas as pd
import numpy as np
from Bio.PDB import *

pdb_name = argv[1]

pdbl = PDBList()
pdbl.retrieve_pdb_file(pdb_name, pdir="path_to_folder/pdb_original/downloaded_pdbs/")

stcrdab_info = pd.read_csv("path_to_folder/pdb_original/stcr_summary_annotated_norepeats.csv", index_col = 0)
alpha = stcrdab_info.loc[pdb_name, "Achain"]
beta = stcrdab_info.loc[pdb_name, "Bchain"]
peptide = stcrdab_info.loc[pdb_name, "antigen_chain"]
mhc1 = stcrdab_info.loc[pdb_name, "mhc_chain1"]
mhc2 = stcrdab_info.loc[pdb_name, "mhc_chain2"]

print(alpha, beta, peptide, mhc1, mhc2)

class chain_select(Select): #function to create TCR + peptide PDB
    def accept_chain(self, chain, alpha = alpha, beta = beta, peptide = peptide, mhc1=mhc1, mhc2=mhc2):
        if chain.full_id[2] == str(alpha) or chain.full_id[2] == str(beta) or chain.full_id[2] == str(peptide) or chain.full_id[2] == str(mhc1) or chain.full_id[2] == str(mhc2):
            return 1
        else:
            return 0

Parser = MMCIFParser()
structure = Parser.get_structure(pdb_name, "path_to_folder/pdb_original/downloaded_pdbs/" + pdb_name + ".cif")

io = PDBIO()
io.set_structure(structure)
io.save("path_to_folder/pdb_original/clean_pdbs/" + pdb_name +".pdb", chain_select())

with open("path_to_folder/pdb_original/clean_pdbs/" + pdb_name +".pdb") as f:
    pdb_file = f.readlines()
    pdb_file = [x.strip() for x in pdb_file]

final_pdb = []

for line in pdb_file:
    if line[0:4] == "ATOM":
        line = list(line)
        if line[21] == alpha:
            new = "A"
        elif line[21] == beta:
            new = "B"
        elif line[21] == peptide:
            new = "P"
        elif line[21] == mhc1:
            new = "M"
        elif line[21] == mhc2:
            new = "N"
        line[21] = new
        new_line = "".join(line)
    else:
        continue
    final_pdb.append(new_line)

np.savetxt("path_to_folder/pdb_original/all_pdbs/" + pdb_name + "_TCR-pMHC.pdb", final_pdb, fmt='%s')
