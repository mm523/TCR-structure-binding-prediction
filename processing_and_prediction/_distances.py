import numpy as np
import pandas as pd
from Bio.PDB import *
import os
from Bio import SeqUtils
from sys import argv
import structure_preprocess_functions as my_preproc

def calculate_distances(arg1, arg2, structure):
    res_numbers = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 112)) + ["111A", "111B",
                    "111C", "111D", "111E", "111F", "112E", "112F", "112D", "112C", "112B", "112A", "112"] + list(range(113, 118))
    res_numbers = [str(x) for x in res_numbers]
    p_nums = [str(x) for x in list(range(1,21))]
    CDR1_start = 27
    CDR1_end = 38
    CDR2_start = 56
    CDR2_end = 65
    CDR3_start = 105
    CDR3_end = 117

    structure_id = arg1.split("_TCR-pMHC.pdb")[0].split(".pdb")[0]
    filename = "renumbered_pdb/" + structure_id + "_renumbered.pdb"
    folder = arg2
    if not os.path.exists(folder+"Results_by_pdb/distances/"):
        os.mkdir(folder+"Results_by_pdb/distances/")

    pep_chain = "P"
    alpha_chain = "A"
    beta_chain = "B"

    ##make PDB with only peptide and TCR to localise the PyRosetta search

    all_atom_distances_alpha, match_table_alpha = my_preproc.calculate_pairwise_distances(alpha_chain, "alpha", structure)
    all_atom_distances_beta, match_table_beta = my_preproc.calculate_pairwise_distances(beta_chain, "beta", structure)

    all_atom_distances_alpha.to_csv(folder+"Results_by_pdb/distances/" + structure_id + "_all_distances_alpha.csv")
    match_table_alpha.to_csv(folder+"Results_by_pdb/distances/" + structure_id + "_distance_table_alpha.csv")
    all_atom_distances_beta.to_csv(folder+"Results_by_pdb/distances/" + structure_id + "_all_distances_beta.csv")
    match_table_beta.to_csv(folder+"Results_by_pdb/distances/" + structure_id + "_distance_table_beta.csv")

    col_names = []
    for row in range(1, 21):
        row = str(row)
        for col in res_numbers:
            col = str(col)
            col_names.append(str(row + "-" + col))

    alpha_unravel = np.array(match_table_alpha).reshape(-1)
    beta_unravel = np.array(match_table_beta).reshape(-1)
    alpha_col_names = [str(x + "_alpha") for x in col_names]
    beta_col_names = [str(x + "_beta") for x in col_names]
#    print(alpha_unravel)
#    print(alpha_col_names)
    alpha_unravel = pd.DataFrame(alpha_unravel.reshape(-1, len(alpha_unravel)), columns = alpha_col_names)
    beta_unravel = pd.DataFrame(beta_unravel.reshape(-1, len(beta_unravel)), columns = beta_col_names)
    linear_distances = pd.concat([alpha_unravel, beta_unravel], axis = 1) #unravels horizontally (each row next to the other)
    linear_distances.to_csv(folder + "outputs_for_model/" + structure_id + "_linearised_distances.csv")

    return linear_distances
