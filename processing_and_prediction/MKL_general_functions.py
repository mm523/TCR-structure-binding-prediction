import pandas as pd
import numpy as np
from MKLpy import metrics
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn import preprocessing, model_selection, decomposition
from MKLpy.model_selection import cross_val_score
from MKLpy.algorithms import EasyMKL
from MKLpy.generators import RBF_generator, HPK_generator
from sklearn.svm import SVC
from itertools import product
import MKL_my_functions1 as my_mkl
import MKL_general_functions as gen_fn
import sklearn.metrics as skmetrics
import matplotlib.pyplot as plt
import argparse
import os
import time
import random


def find_interesting_gammas(X_train, Y_train):
    ### applying heuristics to find best gammas to use in RBF
    # basically you find distances between positives and negatives and then quantiles
    X_train_pos = X_train[Y_train == 1, :]
    X_pos_selec = random.choices(list(range(0, np.shape(X_train_pos)[0])), k = min(np.shape(X_train_pos)[0], 1000))
    X_pos_subset = X_train_pos[X_pos_selec,:]
    X_train_neg = X_train[Y_train == 0, :]
    X_neg_selec = random.choices(list(range(0, np.shape(X_train_neg)[0])), k = min(np.shape(X_train_neg)[0], 1000))
    X_neg_subset = X_train_neg[X_neg_selec,:]
    distances = []
    for i in range(0, np.shape(X_pos_subset)[0]):
        for j in range(0, np.shape(X_neg_subset)[0]):
            dist = np.linalg.norm(X_pos_subset[i,:]-X_neg_subset[j, :])
            distances.append(dist)
    sigmas = []
    gammas= []
    # print(distances)
    for perc in [1, 2, 5, 50, 95, 98, 99]:
        sigma = np.percentile(distances, perc)
        sigmas.append(sigma)
        gamma = 1/(2*sigma**2)
        gammas.append(gamma)
    return gammas

def find_interesting_gammas_no0(X_train, Y_train):
    ### applying heuristics to find best gammas to use in RBF
    # basically you find distances between positives and negatives and then quantiles
    X_train_pos = X_train[Y_train == 1, :]
    X_pos_selec = random.choices(list(range(0, np.shape(X_train_pos)[0])), k = min(np.shape(X_train_pos)[0], 1000))
    X_pos_subset = X_train_pos[X_pos_selec,:]
    X_train_neg = X_train[Y_train == 0, :]
    X_neg_selec = random.choices(list(range(0, np.shape(X_train_neg)[0])), k = min(np.shape(X_train_neg)[0], 1000))
    X_neg_subset = X_train_neg[X_neg_selec,:]
    distances = []
    for i in range(0, np.shape(X_pos_subset)[0]):
        for j in range(0, np.shape(X_neg_subset)[0]):
            dist = np.linalg.norm(X_pos_subset[i,:]-X_neg_subset[j, :])
            distances.append(dist)
    sigmas = []
    gammas= []
    # print(distances)
    for perc in [1, 2, 5, 50, 95, 98, 99]:
        sigma = np.percentile(distances, perc)
        if sigma == 0:
            continue
        else:
            sigmas.append(sigma)
            gamma = 1/(2*sigma**2)
            gammas.append(gamma)
    return gammas

def prepare_dataset(dataset, pick, good_templates):
    # print(good_templates)
    label_list = []
    to_drop = []
    dataset = dataset.sort_values(by = ["structure_id"]).reset_index(drop=True)
    dataset = dataset[dataset.structure_id != "decoy158_decoy_original"]
    dataset = dataset[dataset.structure_id != "decoy158_2_decoy_original"]
    dataset = dataset[dataset.structure_id != "decoy158_decoy_original_2"]
    for i in range(0, np.shape(dataset)[0]):
        line = dataset.iloc[i]
        # print(line)
        # print(line.structure_id)
        # print(line.structure_id in good_templates)
        if (line.binding == 0) and (line.structure_id in good_templates):
            label_list.append("non-binder")
        elif (line.structure_id in good_templates) and (line.binding == 1):
            label_list.append("binder")
        else:
            to_drop.append(dataset.index.values[i])
    # print("dropping")
    # print(to_drop)
    dataset = dataset.drop(to_drop, axis = 0)
    ids = list(dataset.structure_id)
    dataset = dataset.drop("structure_id", axis = 1)
    dataset = dataset.drop("binding", axis = 1)
    size = np.shape(dataset)[0]
    print("dataset size: ", size)
    print("positives: ", sum([int(x == "binder") for x in label_list]))
    print("negatives: ", sum([int(x == "non-binder") for x in label_list]))
    for label in dataset.columns.values:
        if "structure_type" in label:
            dataset = dataset.drop(label, axis = 1)
    if pick != 1:
        dataset, hold_out = model_selection.train_test_split(dataset, train_size = pick, random_state = 0)
        label_list, hold_out_labels = model_selection.train_test_split(label_list, train_size = pick, random_state = 0)

    return dataset, label_list, size, ids

def prepare_datasets_without_var_j(dataset, pick, good_templates, dataset_name):
    label_list = []
    to_drop = []
    dataset = dataset.sort_values(by = ["structure_id"]).reset_index(drop=True)
    dataset.index = list(dataset.structure_id)
    dataset = dataset[dataset.structure_id != "decoy158_decoy_original"]
    dataset = dataset[dataset.structure_id != "decoy158_2_decoy_original"]
    dataset = dataset[dataset.structure_id != "decoy158_decoy_original_2"]
    for i in range(0, np.shape(dataset)[0]):
        line = dataset.iloc[i]
        if (line.binding == 0) and (line.structure_id in good_templates):
            label_list.append("non-binder")
        elif (line.structure_id in good_templates) and (line.binding == 1):
            label_list.append("binder")
        else:
            to_drop.append(dataset.index.values[i])
    dataset = dataset.drop(to_drop, axis = 0)
    dataset = dataset.drop("structure_id", axis = 1)
    col_names = dataset.columns.values
    # if ds_j == dataset_name:
    #     var_j_data = dataset[var_j]
    #     dataset_no_varj = dataset.drop(var_j, axis = 1)
    # else:
    #     dataset_no_varj = dataset
    dataset_half = dataset.sample(frac = 0.5, axis = 1)
    size = np.shape(dataset)[0]
    # print("dataset size: ", size)
    # print("positives: ", sum([int(x == "binder") for x in label_list]))
    # print("negatives: ", sum([int(x == "non-binder") for x in label_list]))
    for label in dataset_half.columns.values:
        if "structure_type" in label:
            dataset_half = dataset_half.drop(label, axis = 1)
    if pick != 1:
        dataset, hold_out = model_selection.train_test_split(dataset, train_size = pick, random_state = 0)
        label_list, hold_out_labels = model_selection.train_test_split(label_list, train_size = pick, random_state = 0)

    # if ds_j == dataset_name:
    #     # print("adding var_j to this dataset")
    #     print(var_j)
    #     dataset_plus_varj = pd.concat([dataset_half, var_j_data], axis = 1)

    return dataset_half, label_list, size #dataset_plus_varj,

def prepare_dataset1(dataset, pick, good_templates): ##also extract structure_id here to take out structures to mutate
    label_list = []
    to_drop = []
    dataset = dataset.sort_values(by = ["structure_id"]).reset_index(drop=True)
    # dataset = dataset[dataset.structure_id != "decoy158_decoy_original"]
    # dataset = dataset[dataset.structure_id != "decoy158_2_decoy_original"]
    # dataset = dataset[dataset.structure_id != "decoy158_decoy_original_2"]
    for i in range(0, np.shape(dataset)[0]):
        line = dataset.iloc[i]
        if "decoy_original" in line.structure_id:
            label_list.append("non-binder")
        elif "model_original" in line.structure_id:
            label_list.append("binder")
        elif "B10x" in line.structure_id:
            if "decoy" in line.structure_id and line.structure_id in good_templates:
                label_list.append("non-binder")
            elif line.structure_id in good_templates:
                label_list.append("binder")
            else:
                to_drop.append(dataset.index.values[i])
        else:
            to_drop.append(dataset.index.values[i])
    dataset = dataset.drop(to_drop, axis = 0)
    structure_ids = list(dataset.structure_id)
    dataset = dataset.drop("structure_id", axis = 1)
    size = np.shape(dataset)[0]
    print("dataset size: ", size)
    print("positives: ", sum([int(x == "binder") for x in label_list]))
    print("negatives: ", sum([int(x == "non-binder") for x in label_list]))
    for label in dataset.columns.values:
        if "structure_type" in label:
            dataset = dataset.drop(label, axis = 1)
    if pick != 1:
        dataset, hold_out = model_selection.train_test_split(dataset, train_size = pick, random_state = 0)
        label_list, hold_out_labels = model_selection.train_test_split(label_list, train_size = pick, random_state = 0)

    return dataset, label_list, size, structure_ids

def kernel_generators(X_train, Y_train, X_test, kernel_type):
    if kernel_type == "rbf":
        gammas = find_interesting_gammas(X_train, Y_train)
        Ktr = RBF_generator(X_train, gamma = gammas, cache=True)
        Kte = RBF_generator(X_test, X_train, gammas, cache=True)
    elif kernel_type == "hpk":
        Ktr = HPK_generator(X_train, degrees=range(1,30), cache=True)
        Kte = HPK_generator(X_test, X_train, degrees=range(1,30), cache=True)
    else:
        print("This kernel: " + kernel_type + "is not yet supported by this script")
    return Ktr, Kte, gammas

def kernel_generators_train(X_train, Y_train, kernel_type):
    if kernel_type == "rbf":
        gammas = find_interesting_gammas(X_train, Y_train)
        Ktr = RBF_generator(X_train, gamma = gammas, cache=True)
    elif kernel_type == "hpk":
        Ktr = HPK_generator(X_train, degrees=range(1,30), cache=True)
    else:
        print("This kernel: " + kernel_type + "is not yet supported by this script")
    return Ktr, gammas

def kernel_generators_train_no0(X_train, Y_train, kernel_type):
    if kernel_type == "rbf":
        gammas = find_interesting_gammas_no0(X_train, Y_train)
        Ktr = RBF_generator(X_train, gamma = gammas, cache=True)
    elif kernel_type == "hpk":
        Ktr = HPK_generator(X_train, degrees=range(1,30), cache=True)
    else:
        print("This kernel: " + kernel_type + "is not yet supported by this script")
    return Ktr, gammas

def kernel_generators_test(X_test, X_train, kernel_type, gammas):
    if kernel_type == "rbf":
        Kte = RBF_generator(X_test, X_train, gamma = gammas, cache=True)
    elif kernel_type == "hpk":
        Kte = HPK_generator(X_test, X_train, degrees=range(1,30), cache=True)
    else:
        print("This kernel: " + kernel_type + "is not yet supported by this script")
    return Kte

def kernel_preparation(dataset, pca_status, external_train_index, external_test_index, Y_labels):
    dataset = dataset.loc[:, (dataset != dataset.iloc[0]).any()]
    dataset = np.array(dataset)
    if external_train_index and external_test_index:
        X_train = dataset[external_train_index]
        print("X_train size: ", np.shape(X_train))
        X_test = dataset[external_test_index]
        Y_train = np.array(Y_labels)[external_train_index]
        Y_test = np.array(Y_labels)[external_test_index]
        X_train = imp.fit_transform(X_train)
        X_train = scaler.fit_transform(X_train)
        X_test = imp.transform(X_test)
        X_test = scaler.transform(X_test)
        if pca_status == "yes":
            X_train = pca.fit_transform(X_train)
            X_test = pca.transform(X_test)
        else:
            pca = None
        X_train = preprocessing.normalize(X_train)
        X_test = preprocessing.normalize(X_test)
    return X_train, Y_train, X_test, Y_test, imp, scaler, pca

def kernel_preparation_train(dataset, pca_status, Y_labels, imp, scaler, pca):
    X_train = dataset
    print("X_train size: ", np.shape(X_train))
    Y_train = np.array(Y_labels)
    X_train = imp.fit_transform(X_train)
    X_train = scaler.fit_transform(X_train)
    if pca_status == "yes":
        X_train = pca.fit_transform(X_train)
    else:
        pca = None
    X_train = preprocessing.normalize(X_train)
    X_test = []
    Y_test = []
    return X_train, Y_train, imp, scaler, pca

def kernel_preparation_test(dataset, pca_status, fit_imputer, fit_scaler, fit_pca = None):
    X_test = dataset
    X_test = fit_imputer.transform(X_test)
    X_test = fit_scaler.transform(X_test)
    if pca_status == "yes":
        X_train = fit_pca.transform(X_train)
    X_test = preprocessing.normalize(X_test)
    return X_test

def find_optimal_params(K_combo_tr, Y_train, cv_int, score_type, lam_values, C_values):
    CV_results = {"lam" : [], "C" : [], "scores" : [], "threshold" : []}

    for lam, C in product(lam_values, C_values):
        print(lam, C)
        svm = SVC(C=C)
        mkl = EasyMKL(lam=lam, learner=svm)
        scores, thresholds,precisions, recalls = my_mkl.my_cross_val_score(K_combo_tr, Y_train, mkl, n_folds=cv_int, scoring=score_type)
        CV_results["lam"].append(lam)
        CV_results["C"].append(C)
        CV_results["scores"].append(np.mean(scores))
        CV_results["threshold"].append(np.mean(scores))
        # except:
        #     print(lam, C)
        #     CV_results["lam"].append(lam)
        #     CV_results["C"].append(C)
        #     CV_results["scores"].append(np.nan)
        #     CV_results["threshold"].append(np.nan)
    CV_results = pd.DataFrame.from_dict(CV_results)
    print(CV_results)
    chosen_lam, chosen_C, chosen_t = CV_results.iloc[CV_results["scores"].idxmax()]["lam"], CV_results.iloc[CV_results["scores"].idxmax()]["C"], CV_results.iloc[CV_results["scores"].idxmax()]["threshold"]
    return chosen_lam, chosen_C, chosen_t
