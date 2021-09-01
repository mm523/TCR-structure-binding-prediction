import pickle
import pandas as pd
import numpy as np
from MKLpy import metrics
from sklearn.impute import SimpleImputer
from MKLpy.model_selection import cross_val_score
from MKLpy.algorithms import EasyMKL
from MKLpy.generators import RBF_generator, HPK_generator
from sklearn.svm import SVC
from itertools import product
import MKL_my_functions1 as my_mkl
import MKL_general_functions as gen_fn
import sklearn.metrics as skmetrics
import matplotlib.pyplot as plt
import os

###initialise model characteristics
pca_status = "no"
kernel_type = "rbf"
today = "2021-05-27"

def predict_from_datasets(atchley, distances, atr, combo):
    datasets = {"distances" : distances, "atchley":atchley, "fa_atr" : atr}
    dataset_combinations = {"distances" : ["distances"], "atchley" : ["atchley"],
                            "atr" : ["fa_atr"],
                            "atchley-dist" : ["atchley", "distances"],
                            "atchley-atr" : ["atchley", "fa_atr"],
                            "atchley-dist-atr" : ["atchley", "distances", "fa_atr"],
                            "dist-atr" : ["distances", "fa_atr"]}
    my_datasets = dataset_combinations[combo]
    Kernels_te = {}

    for key in my_datasets:
        print("processing single dataset: " + key)
        dataset_results = {}
        dataset = pd.DataFrame(datasets[key])
        imp = pickle.load(open("path_to_folder/trained_models/"+ key + "_" + combo + "_imputer_all-templ_" + today + ".pkl", "rb"))
        scaler = pickle.load((open("path_to_folder/trained_models/" + key + "_" + combo + "_scaler_all-templ_" + today + ".pkl", "rb")))
        with open("path_to_folder/trained_models/rbf_gammas_" + key + "_" + combo + "_all-templ_" + today + ".txt") as f:
            gammas = f.readlines()
        gammas = [float(x) for x in gammas]
        X_train = np.genfromtxt("path_to_folder/trained_models/X_train_" + key + "_" + combo + "_all-templ.csv", delimiter=",")
        X_test = gen_fn.kernel_preparation_test(dataset, pca_status, imp, scaler)
        Kte = gen_fn.kernel_generators_test(X_test, X_train, kernel_type, gammas)
        Kernels_te[key]=Kte

    filename = "path_to_folder/trained_models/PDB_all_templ_" + combo + "_" + today + ".sav"
    trained_model = pickle.load(open(filename, "rb"))

    K_combo_te = []
    for dataset_name in my_datasets:
        K_combo_te = K_combo_te + list(Kernels_te[dataset_name])

    prediction = trained_model.predict(K_combo_te)
    dec_fun = trained_model.decision_function(K_combo_te)

    return prediction, dec_fun
