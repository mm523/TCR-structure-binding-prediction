## this script collect all the info from the results log files and collates them in a single file

import pandas as pd
from sys import argv
import numpy as np
import os
from sys import argv
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')
set = argv[1]

def parse_results(folder, out_name, set):
    results = pd.DataFrame()
    idx1 = 0

    for idx, f in enumerate(os.listdir(folder)):
        id = f.split("_log")[0]
        my_file = id + "_log-good_templ_" + set + ".txt"
        print(my_file)
        idx1 += 1
        template = ""
        dec_function = ""
        classification = ""
        if idx1 == 0 or idx % 1000 == 0:
            print(idx1)
        structure_name = id
        if structure_name != "":
            with open(folder + my_file) as f:
                file = f.readlines()
            for line in file:
                if line.split(":")[0] == "the model says":
                    classification = line.split(":")[1].strip().strip("[]")
                elif line.split(":")[0] == "decision function score":
                    dec_function = line.split(":")[1].strip().strip("[]")
            results.loc[structure_name, "decision_function"] = dec_function
            results.loc[structure_name, "model_answer"] = classification
    results.to_csv(out_name)

b10x_folder = "path_to_folder/B10x/"
dash_folder = "path_to_folder/Dash/"
expt_folder = "path_to_folder/expt/"
atlas_folder = "path_to_folder/atlas/"
# os_folder = "path_to_folder/OS/"
# zhang_folder = "path_to_folder/zhang/"
newVdj_folder = "path_to_folder/newVdj/"
logs_folders = "many_models_logs/"

out_b10x = "path_to_folder/B10x/B10x_prediction_results_gt_" + set + ".csv"
out_dash = "path_to_folder/Dash/Dash_prediction_results_gt_" + set + ".csv"
out_expt = "path_to_folder/expt/expt_prediction_results_gt_" + set + ".csv"
out_atlas = "path_to_folder/atlas/atlas_prediction_results_gt_" + set + ".csv"
# out_os = "path_to_folder/OS/OS_prediction_results_gt_" + set + ".csv"
# out_zhang = "path_to_folder/zhang/zhang_prediction_results_gt_" + set + ".csv"
out_newVdj = "path_to_folder/newVdj/newVdj_prediction_results_gt_" + set + ".csv"

parse_results(b10x_folder + logs_folders, out_b10x, set)
parse_results(dash_folder + logs_folders, out_dash, set)
parse_results(expt_folder + logs_folders, out_expt, set)
parse_results(atlas_folder + logs_folders, out_atlas, set)
# parse_results(os_folder + logs_folders, out_os, set)
# parse_results(zhang_folder + logs_folders, out_zhang, set)
parse_results(newVdj_folder + logs_folders, out_newVdj, set)
