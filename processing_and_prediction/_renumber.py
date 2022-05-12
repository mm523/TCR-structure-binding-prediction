import numpy as np
import pandas as pd
from Bio.PDB import *
import os
from anarci import anarci
from Bio import SeqUtils
from difflib import SequenceMatcher
from sys import argv
import structure_preprocess_functions as my_preproc

def renumbering(arg1, arg2):
    seq_match = SequenceMatcher()
    structure_id = arg1.split("_TCR-pMHC.pdb")[0].split(".pdb")[0]
    filename = "all_pdbs/" + arg1
    folder = arg2

    pep_chain = "P"
    alpha_chain = "A"
    beta_chain = "B"

    Parser = PDBParser()
    structure = Parser.get_structure(structure_id, folder + filename)

    #extract sequences of my structures
    sequence_A, seq_num_A, seq_res_A = my_preproc.extract_sequence(alpha_chain, structure)
    sequence_B, seq_num_B, seq_res_B = my_preproc.extract_sequence(beta_chain, structure)
    sequence_P, seq_num_P, seq_res_P = my_preproc.extract_sequence(pep_chain, structure)

    IMGT_A = anarci([("seq",sequence_A)], "IMGT") #run anarci for renumbering
    IMGT_B = anarci([("seq",sequence_B)], "IMGT")
    IMGT_dict_A = IMGT_A[0][0][0][0] #this is how I get the info I want
    IMGT_dict_B = IMGT_B[0][0][0][0]
    print(IMGT_dict_B)

    alpha_IMGT = my_preproc.IMGT_vs_PDB_numbering(IMGT_dict_A, seq_res_A, seq_num_A)
    beta_IMGT = my_preproc.IMGT_vs_PDB_numbering(IMGT_dict_B, seq_res_B, seq_num_B)

    np.savetxt(folder+"anarci/" + structure_id + "_IMGT_numbering_alpha.csv", alpha_IMGT, delimiter=",", fmt='%s', header="IMGT_num, Res, pdb_num")
    np.savetxt(folder+"anarci/" + structure_id + "_IMGT_numbering_beta.csv", beta_IMGT, delimiter=",", fmt='%s', header="IMGT_num, Res, pdb_num")

    ##import PDB file as text files

    final_pdb = my_preproc.standard_renumbering(sequence_P, alpha_IMGT, beta_IMGT, folder, filename, seq_num_P, pep_chain, alpha_chain, beta_chain)
    np.savetxt(folder+"renumbered_pdb/" + structure_id + "_renumbered.pdb", final_pdb, fmt='%s')
