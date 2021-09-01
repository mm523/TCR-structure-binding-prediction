import pandas as pd
import numpy as np

info = pd.read_excel("path_to_folder/experimental_constructs/Experimental_constructs.xlsx")
print(info)

peptides = list(info.pep_name)

all_experimental_constructs = pd.DataFrame()

for idx in info.index.values:
    construct_line = info.iloc[idx]
    construct = construct_line.construct
    pep_name_or = construct_line.pep_name
    for pep_name in peptides:
        id = construct + "-" + pep_name
        print(id)
        alpha = construct_line.alpha
        beta = construct_line.beta
        pep = info[info.pep_name == pep_name]["peptide"].iloc[0]
        mhc = info[info.pep_name == pep_name]["mhc"].iloc[0]
        mhc_gene = info[info.pep_name == pep_name]["mhc_gene"].iloc[0]
        all_experimental_constructs.loc[id, "construct"] = construct
        all_experimental_constructs.loc[id, "pep_name"] = pep_name
        all_experimental_constructs.loc[id, "mhc_gene"] = mhc_gene

all_experimental_constructs.to_csv("path_to_folder/experimental_constructs/all_expt_constructs.csv")
