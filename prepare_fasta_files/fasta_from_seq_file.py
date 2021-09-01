import pandas as pd
import numpy as np
from sys import argv

def make_fasta(seqs, out_folder):
    for idx in seqs.index.values:
        alpha = seqs.loc[idx, "alpha"]
        beta = seqs.loc[idx, "beta"]
        peptide = seqs.loc[idx, "peptide"]
        mhc = seqs.loc[idx, "mhc"]
        f = open(out_folder + "fasta_files/" + idx + ".fasta", "w+")
        chain_names = ["alpha", "beta", "mhc", "peptide"]
        chains = [alpha, beta, mhc, peptide]
        for chain_name, chain in zip(chain_names, chains):
            print(chain)
            f.write(">" + chain_name + "\n")
            f.write(chain + "\n")
        f.close()

f = "path_to_folder/"

folder_pdb = f + "PDB_set/"
folder_dash = f + "Dash/"
folder_b10x = f + "B10x/"
folder_expt = f + "expt/"
folder_atlas = f + "atlas/"
folder_zhang = f + "zhang/"
folder_zhang1 = f + "zhang1/"
folder_os = f + "OS/"
folder_newVdj = f + "newVdj/"

# print("PDB set: ...")
# seqs_pdb = pd.read_csv(folder_pdb + "PDB_set_all_sequences.csv", index_col = 0).dropna(how = "all", axis = 0).drop_duplicates()
# print("dash: ...")
# seqs_dash = pd.read_csv(folder_dash + "Dash_all_sequences_with_binding.csv", index_col = 0).dropna(how = "all", axis = 0).drop_duplicates()
# print("b10x: ...")
# seqs_b10x = pd.read_csv(folder_b10x + "10x_sequences_with_binding.csv", index_col = 0).dropna(how = "all", axis = 0).drop_duplicates()
# print("expt: ...")
# seqs_expt = pd.read_csv(folder_expt + "expt_constructs_sequences.csv", index_col = 0).dropna(how = "all", axis = 0).drop_duplicates()
# print("atlas: ...")
# seqs_atlas = pd.read_csv(folder_atlas + "atlas_all_sequences.csv", index_col = 0).dropna(how = "all", axis = 0).drop_duplicates()
# print("zhang: ...")
# seqs_zhang = pd.read_csv(folder_zhang + "zhang_all_sequences.csv", index_col = 0).dropna(how = "any", axis = 0).drop_duplicates()
# print("OS: ...")
# seqs_os = pd.read_csv(folder_os + "OS_all_sequences_with_binding.csv", index_col = 0).dropna(how = "any", axis = 0).drop_duplicates()
# print("zhang1: ...")
# seqs_zhang1 = pd.read_csv(folder_zhang1 + "zhang_all_sequences1.csv", index_col = 0).dropna(how = "any", axis = 0).drop_duplicates()
print("newVdj: ...")
seqs_newVdj = pd.read_csv(folder_newVdj + "newVdj_all_sequences.csv", index_col = 0).dropna(how = "any", axis = 0).drop_duplicates()


# make_fasta(seqs_pdb, folder_pdb)
# make_fasta(seqs_dash, folder_dash)
# make_fasta(seqs_b10x, folder_b10x)
# make_fasta(seqs_expt, folder_expt)
# make_fasta(seqs_atlas, folder_atlas)
# make_fasta(seqs_zhang, folder_zhang)
# make_fasta(seqs_os, folder_os)
# make_fasta(seqs_zhang1, folder_zhang1)
make_fasta(seqs_newVdj, folder_newVdj)
