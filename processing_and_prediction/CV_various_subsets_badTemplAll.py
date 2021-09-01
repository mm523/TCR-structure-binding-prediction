### this script performs CV on PDB set to evaluate performance with this latest implementation.
### then, it trains a final model on 70% of the set and tests on 30% of the set (which includes 2nx5, on which I plan on doing the mutation work)
### I will save the model on 70% to try and use for feature selection

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
from sys import argv
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')

"""
About scaling:
        'many elements used in the objective function of a learning algorithm
        (such as the RBF kernel of Support Vector Machines or the l1 and l2
        regularizers of linear models) assume that all features are centered
        around zero and have variance in the same order.'

About arguments:
        I use the arguments that I have established to work best by CV.

"""
##set arguments - these have been validated by CV

pick = 1
combo_to_run = str(argv[1]).strip()
templ_thresh = str(argv[2]).strip()
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

distances = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_distances.csv", index_col = 0)
fa_atr = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_atr.csv", index_col = 0)
fa_elec = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_elec.csv", index_col = 0)
fa_sol = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_sol.csv", index_col = 0)
fa_rep = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_rep.csv", index_col = 0)
atchley = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_AF.csv", index_col = 0)
templ = pd.read_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/PDB_set_badTemplAll_t" + templ_thresh + "_linearised_templates.csv", index_col = 0)

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
good_templates = list(distances.structure_id) #all structures

for idx, key in enumerate(datasets_of_interest): #iterating over the datasets for which I want to create a kernel
    print("###" + key + "####")
    dataset0 = datasets0[key]
    print(dataset0)
    dataset0, label_list0, size0, ids = gen_fn.prepare_dataset(dataset0, 1, good_templates) #setting pick = 1 because I want all PDBs
    print(dataset0)
    dataset = dataset0
    label_list = label_list0
    labels[key] = [int(x == "binder") for x in label_list]
    datasets[key] = np.array(dataset)
    size = np.shape(dataset)[0]

indeces = list(range(0, size))
Y_labels = labels[key] #they are now all the same because we've sorted the datasets first

lam_values = [0]#np.linspace(0, 1, 20) #I could validate for lambda, but it does not seem to be necessary from the Lauriola paper
C_values   = np.logspace(-5, 2, 25)

CV_results = pd.DataFrame()
Kernels_tr = {}
Kernels_te = {}

for idx, (external_train_index, external_test_index) in enumerate(external_fold.split(indeces, indeces)):
    test_ids = np.array(ids)[external_test_index]
    for key in datasets_of_interest: #create kernels for each dataset separately
        print("single dataset results: " + key)
        dataset_results = {}
        dataset = np.array(datasets[key])
        dataset_train = dataset[external_train_index, :]
        dataset_test = dataset[external_test_index, :]
        Y_train = np.array(Y_labels)[external_train_index]
        Y_test = np.array(Y_labels)[external_test_index]
        # print(np.shape(dataset_train), np.shape(Y_train))
        # print(dataset_train)
        # print(Y_train)
        X_train, Y_train, imp, scaler, pca = gen_fn.kernel_preparation_train(dataset_train, pca_status, Y_train, imp, scaler, pca)
        Ktr, gammas = gen_fn.kernel_generators_train(X_train, Y_train, kernel_type)
        print(gammas)
        X_test = gen_fn.kernel_preparation_test(dataset_test, pca_status, imp, scaler)
        Kte = gen_fn.kernel_generators_test(X_test, X_train, kernel_type, gammas)
        Kernels_tr[key]=Ktr
        Kernels_te[key]=Kte


    combo = dataset_combinations[combo_to_run]
    print("Multiple datasets combos - ", combo_to_run, "\n", "datasets included: ", combo)
    K_combo_tr = []
    K_combo_te = []
    for dataset_name in combo:
        K_combo_tr = K_combo_tr + list(Kernels_tr[dataset_name])
        K_combo_te = K_combo_te + list(Kernels_te[dataset_name])

    #internal CV to find optimal model parameters
    chosen_lam, chosen_C, chosen_t = gen_fn.find_optimal_params(K_combo_tr, Y_train, cv_int, score_type, lam_values, C_values)

    svm_best = SVC(C=chosen_C, probability = True)
    mkl_best = EasyMKL(lam = chosen_lam, learner = svm_best)
    mkl_best.fit(K_combo_tr, Y_train) #fit on whole train set
    y_tr_pred = mkl_best.predict(K_combo_tr)
    accuracy_tr, recall_tr, precision_tr, roc_auc_tr = skmetrics.accuracy_score(Y_train, y_tr_pred), skmetrics.recall_score(Y_train, y_tr_pred), skmetrics.precision_score(Y_train, y_tr_pred), skmetrics.roc_auc_score(Y_train, y_tr_pred)
    print ('Accuracy score train: %.4f, recall score train: %.4f, precision score train: %.4f, AUC train: %.4f' % (accuracy_tr, recall_tr, precision_tr, roc_auc_tr))

    y_te_pred = mkl_best.predict(K_combo_te)
    y_te_dec_fun = mkl_best.decision_function(K_combo_te)
    accuracy_te, recall_te, precision_te, roc_auc_te = skmetrics.accuracy_score(Y_test, y_te_pred), skmetrics.recall_score(Y_test, y_te_pred), skmetrics.precision_score(Y_test, y_te_pred), skmetrics.roc_auc_score(Y_test, y_te_pred)
    print ('Accuracy score test: %.4f, recall score test: %.4f, precision score test: %.4f, AUC test: %.4f' % (accuracy_te, recall_te, precision_te, roc_auc_te))

    CV_results.loc[idx, "accuracy_train"] = accuracy_tr
    CV_results.loc[idx, "precision_train"] = recall_tr
    CV_results.loc[idx, "recall_train"] = recall_tr
    CV_results.loc[idx, "accuracy_test"] = accuracy_te
    CV_results.loc[idx, "precision_test"] = precision_te
    CV_results.loc[idx, "recall_test"] = recall_te
    CV_results.loc[idx, "auc_test"] = roc_auc_te
    CV_results.loc[idx, "AP_test"] = skmetrics.average_precision_score(Y_test, mkl_best.decision_function(K_combo_te))
    results_by_structure = pd.DataFrame([y_te_dec_fun, y_te_pred, Y_test]).transpose()
    results_by_structure.columns = ["dec_fun", "prediction", "real_answer"]
    results_by_structure.index = test_ids
    results_by_structure.to_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/CV/PDB_CV_results-" + str(combo_to_run) + "_t" + templ_thresh + "-roc-singlePred-fold" + str(idx) + "_" + today + ".csv")
    # print(results_by_structure)

print(CV_results)
CV_results.to_csv("path_to_folder/PDB_set_badTemplAll_t" + templ_thresh + "/CV/PDB_CV_results-" + str(combo_to_run) + "_t" + templ_thresh + "-roc_" + today + ".csv")
