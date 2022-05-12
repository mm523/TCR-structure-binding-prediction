import pandas as pd
import numpy as np
from Bio.PDB import *
import os
from itertools import product
from sys import argv
import structure_preprocess_functions as my_preproc

def atchley_factors(arg1, arg2):
    structure_id = arg1.split("_TCR-pMHC.pdb")[0].split(".pdb")[0]
    filename = "renumbered_pdb/" + structure_id + "_renumbered.pdb"
    folder = "/home/regmili/Scratch/" + arg2
    folder1 = folder+ "Results_by_pdb/sequences/"
    if not os.path.exists(folder+"Results_by_pdb/sequences/"):
        os.mkdir(folder+"Results_by_pdb/sequences/")

    aa_letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    Atchley_factors = pd.read_csv("path_to_folder/Atchley_factors_table.csv", encoding='utf-8-sig')
    chains = {}
    chains["peptide"] = "P"
    chains["alpha"] = "A"
    chains["beta"] = "B"
    Parser = PDBParser()
    structure = Parser.get_structure(structure_id, folder + filename)

    sequence_results = my_preproc.extract_CDR_and_peptide_sequences(structure_id, structure, chains)
    #sequence_results["structure_id"] = structure_id

    AFs = ["AF1", "AF2", "AF3", "AF4", "AF5"]
    column_names = list(sequence_results.columns.values)
    new_column_names_AF = []
    for name, AF in product(column_names, AFs):
         new_column_names_AF.append(name + "_" + str(AF))

    sequence_results_vector_encoding = pd.DataFrame(columns = new_column_names_AF)

    for j in range(0, np.shape(sequence_results)[1]):
        col_name = sequence_results.columns.values[j].strip()
        aa = sequence_results.loc[structure_id, col_name]
        if pd.isna(aa) == False:
            for AF in AFs:
                new_value = Atchley_factors.loc[Atchley_factors["aa"] == aa, str(AF)].values[0]
    #            print(new_value)
                sequence_results_vector_encoding.loc[structure_id,col_name + "_" + str(AF)] = float(new_value.replace('\u2212', '-'))

    sequence_results_vector_encoding.to_csv(folder1 + structure_id + "_sequence_as_AFs.csv")
    sequence_results_vector_encoding.to_csv(folder + "outputs_for_model/" + structure_id + "_AFs.csv")

    return np.array(sequence_results_vector_encoding)
