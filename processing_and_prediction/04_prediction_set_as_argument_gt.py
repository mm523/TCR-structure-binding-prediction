from sys import argv
import _predict_from_diff_sets_gt as _predict
import pandas as pd

structure_id = argv[1].split("_TCR-pMHC.pdb")[0]
folder = "path_to_folder/" + argv[2]
combo = argv[3]

print("------------------------------------")
print()
print("now predicting whether structure is a binder...")
distances = pd.read_csv(folder + "outputs_for_model/" + structure_id + "_linearised_distances.csv", index_col = 0)
fa_atr = pd.read_csv(folder + "outputs_for_model/" + structure_id + "_linearised_atr.csv", index_col = 0)
# fa_elec = pd.read_csv(folder + "outputs_for_model/" + structure_id + "_linearised_elec.csv", index_col = 0)
atchley = pd.read_csv(folder + "outputs_for_model/" + structure_id + "_AFs.csv", index_col = 0)
templates = pd.read_csv(folder + "outputs_for_model/" + structure_id +  "_template_info.csv", index_col=0)
prediction, dec_fun = _predict.predict_from_datasets(atchley, distances, fa_atr, combo) ## run here the script that predicts given the best existing model
print("the model says: " , prediction)
print("decision function score:", dec_fun)
print()
print()
print("You're welcome :) ")
