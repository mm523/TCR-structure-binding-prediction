from sys import argv
import pandas as pd
import numpy as np
import os

structure_id = argv[1].split("_TCR-pMHC.pdb")[0]
folder = "path_to_folder/" + argv[2]
folder1 = folder + "model_outputs/"

chains = ["pep", "mhc", "alpha", "beta"]
terms = ["avg", "min", "max"]
template = ["complex-templates", "pmhc-templates"]
structure_names = []
results = pd.DataFrame()

for template_type in template:
    template_file = pd.read_csv(folder1 + structure_id + "-" + template_type + ".csv", index_col=0)
#        print(template_file)
    for column_name in template_file.columns.values:
        chain = column_name.split("_")[0]
        data = template_file[column_name]
        results.loc[structure_id, chain + "_avg_" + template_type] = np.mean(data)
        results.loc[structure_id, chain + "_min_" + template_type] = np.min(data)
        results.loc[structure_id, chain + "_max_" + template_type] = np.max(data)

results.to_csv(folder + "outputs_for_model/" + structure_id + "_template_info.csv")
