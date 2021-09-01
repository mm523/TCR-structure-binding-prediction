import pandas as pd
import numpy as np
import os
from sys import argv

s = argv[1]
folder = "path_to_folder/" + s + "/outputs_for_model/"

structure_names = []
for filename in os.listdir(folder):
    if "distances" in filename:
        structure_id = filename.split("_linearised")[0].split("_AF")[0].split("_template")[0]
        structure_names.append(structure_id)
structure_names = set(structure_names)

atr_all = pd.DataFrame()
elec_all = pd.DataFrame()
dist_all = pd.DataFrame()
rep_all = pd.DataFrame()
sol_all = pd.DataFrame()
AF_all = pd.DataFrame()
templ_all = pd.DataFrame()

for structure in structure_names:
    try:
        atr = pd.read_csv(folder + structure + "_linearised_atr.csv", index_col=0)
        elec = pd.read_csv(folder + structure + "_linearised_elec.csv", index_col=0)
        rep = pd.read_csv(folder + structure + "_linearised_rep.csv", index_col=0)
        sol = pd.read_csv(folder + structure + "_linearised_sol.csv", index_col=0)
        dist = pd.read_csv(folder + structure + "_linearised_distances.csv", index_col=0)
        af = pd.read_csv(folder + structure + "_AFs.csv", index_col=0)
    except:
        continue
    try:
        templ = pd.read_csv(folder + structure + "_template_info.csv", index_col = 0)
    except:
        templ = pd.DataFrame()

    af["structure_id"] = structure
    AF_all = pd.concat([AF_all, af])
    atr["structure_id"] = structure
    atr_all = pd.concat([atr_all, atr])
    elec["structure_id"] = structure
    elec_all = pd.concat([elec_all, elec])
    rep["structure_id"] = structure
    rep_all = pd.concat([rep_all, rep])
    sol["structure_id"] = structure
    sol_all = pd.concat([sol_all, sol])
    dist["structure_id"] = structure
    dist_all = pd.concat([dist_all, dist])
    templ["structure_id"] = structure
    templ_all = pd.concat([templ_all, templ])


dist_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_distances.csv")
atr_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_atr.csv")
elec_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_elec.csv")
sol_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_sol.csv")
rep_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_rep.csv")
AF_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_AF.csv")
templ_all.to_csv("path_to_folder/" + s + "/" + s + "_linearised_templates.csv")
