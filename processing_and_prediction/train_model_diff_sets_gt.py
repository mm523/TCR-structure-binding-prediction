## train the model with PDB data
# MKL_my_functions1 contains my edited fucntions from MKLpy
# MKL_general_functions contains the preprocessing steps to make the kernels

import numpy as np
import pandas as pd
from itertools import product
from sklearn.impute import SimpleImputer
from sklearn import preprocessing, decomposition, model_selection
from sklearn.svm import SVC
import sklearn.metrics as skmetrics
from MKLpy import metrics
from MKLpy.model_selection import cross_val_score
from MKLpy.algorithms import EasyMKL
from MKLpy.generators import RBF_generator, HPK_generator
import MKL_my_functions1 as my_mkl
import MKL_general_functions as gen_fn
import os
import pickle
from datetime import datetime
from sys import argv
today = datetime.today().strftime('%Y-%m-%d')

"""
About scaling:
        "many elements used in the objective function of a learning algorithm
        (such as the RBF kernel of Support Vector Machines or the l1 and l2
        regularizers of linear models) assume that all features are centered
        around zero and have variance in the same order."

About arguments:
        I use the arguments that I have established to work best by CV.

"""
##set arguments - these have been validated by CV

pick = 1
combo_to_run = argv[1]
pca_status = "no"
pca_ncomp = None
scaler_type = "minmax"
kernel_type = "rbf"
cv_int = 10
cv_ext = 10
score_type = ["roc_auc"]

imp = SimpleImputer(missing_values = np.nan, strategy = "median")
if scaler_type == "robust":
    scaler = preprocessing.RobustScaler()
elif scaler_type == "standard":
    scaler = preprocessing.StandardScaler()
elif scaler_type == "minmax":
    scaler = preprocessing.MinMaxScaler()
else:
    print("This kernel: " + kernel_type + "is not yet supported by this script")

pca = decomposition.PCA(n_components=pca_ncomp)

##define what datasets the various combinations need to include
dataset_combinations = {"distances" : ["distances"], "templates" : ["templ"], "atchley" : ["atchley"],
                        "atr" : ["fa_atr"], "elec" : ["fa_elec"], "rep" : ["fa_rep"], "sol" : ["fa_sol"],
                        "atr-sol" : ["fa_atr", "fa_sol"],
                        "atr-sol-elec" : ["fa_atr", "fa_sol", "fa_elec"],
                        "atr-sol-elec-rep" : ["fa_atr", "fa_sol", "fa_elec", "fa_rep"],
                        "atchley-dist" : ["atchley", "distances"],
                        "atchley-atr" : ["atchley", "fa_atr"],
                        "atchley-dist-atr" : ["atchley", "distances", "fa_atr"],
                        "atchley-dist-templates" : ["atchley", "distances", "templ"],
                        "atchley-dist-atr-templates": ["atchley", "distances", "fa_atr", "templ"],
                        "atchley-atr-templates": ["atchley", "fa_atr", "templ"],
                        "dist-atr" : ["distances", "fa_atr"],
                        "dist-atr-sol" : ["distances", "fa_atr", "fa_sol"],
                        "dist-atr-sol-elec" : ["distances", "fa_atr", "fa_sol", "fa_elec"],
                        "dist-atr-sol-elec-rep" : ["distances", "fa_atr", "fa_sol", "fa_elec", "fa_rep"],
                        "dist-atr-templates" : ["atchley", "distances", "fa_atr", "templ"]}

datasets_of_interest = dataset_combinations[combo_to_run]

y_scores_0_all=dict((el, []) for el in dataset_combinations)
y_scores_1_all=dict((el, []) for el in dataset_combinations)

distances = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_distances.csv", index_col = 0)
fa_atr = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_atr.csv", index_col = 0)
fa_elec = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_elec.csv", index_col = 0)
fa_sol = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_sol.csv", index_col = 0)
fa_rep = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_rep.csv", index_col = 0)
atchley = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_AF.csv", index_col = 0)
templ = pd.read_csv("path_to_folder/PDB_set/PDB_set_linearised_templates.csv", index_col = 0)

distances["binding"]= [int("decoy" not in x) for x in distances.structure_id]
fa_atr["binding"]= [int("decoy" not in x) for x in fa_atr.structure_id]
fa_elec["binding"]= [int("decoy" not in x) for x in fa_elec.structure_id]
fa_sol["binding"]= [int("decoy" not in x) for x in fa_sol.structure_id]
fa_rep["binding"]= [int("decoy" not in x) for x in fa_rep.structure_id]
atchley["binding"]= [int("decoy" not in x) for x in atchley.structure_id]
templ["binding"]= [int("decoy" not in x) for x in templ.structure_id]

labels = {}
datasets0 = {"distances" : distances, "fa_atr":fa_atr, "fa_elec":fa_elec, "fa_sol":fa_sol, "fa_rep":fa_rep, "atchley":atchley, "templ":templ}
external_fold = model_selection.ShuffleSplit(n_splits = cv_ext, test_size= 0.3)
datasets = {}
fold = -1

good_templates = list(templ[(templ["total_max_complex-templates"] > 60) & (templ["pep_max_pmhc-templates"] > 45)]["structure_id"])

for idx, key in enumerate(datasets_of_interest): #iterating over the datasets for which I want to create a kernel
    print("###" + key + "_" + combo_to_run + "####")
    dataset0 = datasets0[key]
    dataset0, label_list0, size0, ids = gen_fn.prepare_dataset(dataset0, 1, good_templates) #setting pick = 1 because I want all PDBs
    print(dataset0)
    dataset = dataset0
    label_list = label_list0
    labels[key] = [int(x == "binder") for x in label_list]
    datasets[key] = np.array(dataset)
    size = np.shape(dataset)[0]

indeces = list(range(0, np.shape(datasets[key])[0]-1))
Y_labels = labels[key] #they are now all the same because we"ve sorted the datasets first

lam_values = [0]#np.linspace(0, 1, 20) #I could validate for lambda, but it does not seem to be necessary from the Lauriola paper
C_values   = np.logspace(-5, 2, 25)

all_results = pd.DataFrame()
Kernels_tr = {}

for key in datasets_of_interest: #create kernels for each dataset separately
    print("single dataset results: " + key)
    dataset_results = {}
    dataset = pd.DataFrame(datasets[key])
    X_train, Y_train, imp, scaler, pca = gen_fn.kernel_preparation_train(dataset, pca_status, Y_labels, imp, scaler, pca)
    # save train to be used for creation of kernel of test set
    np.savetxt("path_to_folder/trained_models/X_train_" + key + "_" + combo_to_run + "_good-templ.csv", X_train, delimiter=",")
    Ktr, gammas = gen_fn.kernel_generators_train(X_train, Y_train, kernel_type)
    Kernels_tr[key]=Ktr
    #use pickle to save imputer, scaler and gammas
    with open("path_to_folder/trained_models/rbf_gammas_" + key + "_" + combo_to_run + "_good-templ_" + today + ".txt","w") as f:
        f.writelines([str(x) + "\n" for x in gammas])
    with open("path_to_folder/trained_models/" + key + "_" + combo_to_run + "_imputer_good-templ_" + today + ".pkl", "wb") as f1:
        pickle.dump(imp, f1)
    with open("path_to_folder/trained_models/" + key + "_" + combo_to_run + "_scaler_good-templ_" + today + ".pkl", "wb") as f2:
        pickle.dump(scaler, f2)

combo = dataset_combinations[combo_to_run]
print("Multiple datasets combos - ", combo_to_run, "\n", "datasets included: ", combo)
K_combo_tr = []
for dataset_name in combo:
    K_combo_tr = K_combo_tr + list(Kernels_tr[dataset_name])

#internal CV to find optimal model parameters
chosen_lam, chosen_C, chosen_t = gen_fn.find_optimal_params(K_combo_tr, Y_train, cv_int, score_type, lam_values, C_values)

svm_best = SVC(C=chosen_C, probability = True)
mkl_best = EasyMKL(lam = chosen_lam, learner = svm_best)
mkl_best.fit(K_combo_tr, Y_train) #fit on whole train set
print("weights: ", mkl_best.solution.weights)

y_tr_pred = mkl_best.predict(K_combo_tr)
accuracy_tr, recall_tr, precision_tr, roc_auc_tr = skmetrics.accuracy_score(Y_train, y_tr_pred), skmetrics.recall_score(Y_train, y_tr_pred), skmetrics.precision_score(Y_train, y_tr_pred), skmetrics.roc_auc_score(Y_train, y_tr_pred)
print ("Accuracy score train: %.4f, recall score train: %.4f, precision score train: %.4f, AUC train: %.4f" % (accuracy_tr, recall_tr, precision_tr, roc_auc_tr))

with open("path_to_folder/trained_models/PDB_good_templ_" + combo_to_run + "_" + today + ".sav", "wb") as f3:
    pickle.dump(mkl_best, f3) #save trained model
