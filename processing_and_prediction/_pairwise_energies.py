import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
from sys import argv
import structure_preprocess_functions as my_preproc

def energies_of_interest(arg1, arg2):
    structure_id = arg1.split("_TCR-pMHC.pdb")[0].split(".pdb")[0]
    folder = "/home/regmili/Scratch/" + arg2
    folder1 = folder+ "Results_by_pdb/energies/"

    structure_names = []
    i = 0

    alpha_atr, beta_atr, linear_atr = my_preproc.get_energy_table(folder1, structure_id, "fa_atr")
    alpha_elec, beta_elec, linear_elec = my_preproc.get_energy_table(folder1, structure_id, "fa_elec")
    alpha_rep, beta_rep, linear_rep = my_preproc.get_energy_table(folder1, structure_id, "fa_rep")
    alpha_sol, beta_sol, linear_sol = my_preproc.get_energy_table(folder1, structure_id, "fa_sol")

    alpha_atr.to_csv(folder1 +  structure_id + "_energy_table_alpha_fa_atr.csv")
    beta_atr.to_csv(folder1 +  structure_id + "_energy_table_beta_fa_atr.csv")
    alpha_elec.to_csv(folder1 +  structure_id + "_energy_table_alpha_fa_elec.csv")
    beta_elec.to_csv(folder1 +  structure_id + "_energy_table_beta_fa_elec.csv")
    alpha_rep.to_csv(folder1 +  structure_id + "_energy_table_alpha_fa_rep.csv")
    beta_rep.to_csv(folder1 +  structure_id + "_energy_table_beta_fa_rep.csv")
    alpha_sol.to_csv(folder1 +  structure_id + "_energy_table_alpha_fa_sol.csv")
    beta_sol.to_csv(folder1 +  structure_id + "_energy_table_beta_fa_sol.csv")
    linear_atr.to_csv(folder + "outputs_for_model/" + structure_id + "_linearised_atr.csv")
    linear_elec.to_csv(folder + "outputs_for_model/" + structure_id + "_linearised_elec.csv")
    linear_rep.to_csv(folder + "outputs_for_model/" + structure_id + "_linearised_rep.csv")
    linear_sol.to_csv(folder + "outputs_for_model/" + structure_id + "_linearised_sol.csv")

    return linear_atr, linear_elec, linear_rep, linear_sol
