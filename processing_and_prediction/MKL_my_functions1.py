from sklearn.metrics import accuracy_score,roc_auc_score, precision_recall_curve, f1_score, fbeta_score, make_scorer, recall_score, average_precision_score
from sklearn.model_selection import StratifiedKFold as KFold
from sklearn import model_selection as skms
import matplotlib.pyplot as plt
import numpy as np
import torch


def __def_score__(score):
    '''internal check'''
#    score = score.lower()
    if score[0] == 'roc_auc':
    	return make_scorer(roc_auc_score), 'decision_function'
    elif score[0] == 'avg_prec':
    	return make_scorer(average_precision_score), 'decision_function'
    elif score[0] == 'accuracy':
    	return make_scorer(accuracy_score), 'predict'
    elif score[0] == 'f_score':
    	return make_scorer(f1_score), 'predict'
    elif score[0] == 'recall':
    	return make_scorer(recall_score), 'predict'
    elif score[0] == 'fbeta':
        beta = score[1]
        return make_scorer(fbeta_score, beta=beta, average="binary"), 'predict'
    else:
    	raise ValueError('%s is not a valid metric. Valid metrics are \'roc_auc\', \'accuracy\', or \'f_score\'.' % score)

def plot_precision_recall_vs_threshold(precisions, recalls, thresholds):
    """
    Modified from:
    Hands-On Machine learning with Scikit-Learn
    and TensorFlow; p.89
    """
    plt.figure(figsize=(8, 8))
    plt.title("Precision and Recall Scores as a function of the decision threshold")
    plt.plot(thresholds, precisions[:-1], "b--", label="Precision")
    plt.plot(thresholds, recalls[:-1], "g-", label="Recall")
    plt.ylabel("Score")
    plt.xlabel("Decision Threshold")
    plt.legend(loc='best')

def my_cross_val_score(KL, Y, estimator, cv=None, n_folds=3, scoring='accuracy', random_state=None, shuffle=True):
    '''performs the cross validation'''
    if scoring[0] == "fbeta":
        beta = scoring[1]
        scorer = make_scorer(fbeta_score, beta=beta)
        f = 'predict'
    else:
        scorer, f = __def_score__(scoring)
#    print(scorer)
    f = getattr(estimator,f)
    n = len(Y)
    cv   = cv or KFold(n_folds, random_state=random_state, shuffle=shuffle)
    results = []
    thresholds = []
    precisions = []
    recalls = []
    for train,test in cv.split(Y,Y):
        KLtr = [K[train][:,train] for K in KL]
        KLte = [K[test ][:,train] for K in KL]
        clf = estimator.fit(KLtr,Y[train])
        y_scores = clf.decision_function(KLte)
        p, r, threshold = precision_recall_curve(Y[test], y_scores)
        p = p[:len(p)-1]
        r = r[:len(r)-1]
#        print("---")
#        for x in p,r,threshold:
#            print(np.shape(x))
#        print(r)
#        print("max recall: ", np.max(r))
        max_r = np.where(r > 0.75)
#        print(max_r)
        p_max_r = p[max_r]
#        print(r[max_r], p[max_r], threshold[max_r])
#        print(np.shape(p), np.shape(p_max_r))
        max_p = np.where(p == np.max(p_max_r))
#        print("max precision: ", np.max(p_max_r))
#        print(max_p)
        best_recall_threshold = threshold[max_p]
#        print(best_recall_threshold)

#        print(p, r, thresholds)
        score = scorer(clf, KLte, Y[test])
        results.append(score)
        thresholds.append(best_recall_threshold[0])
        precisions.append(p[max_p][0])
        recalls.append(r[max_p][0])
        # print(results)
        # print(thresholds)
    return results, thresholds, precisions, recalls





def train_test_split(KL, Y, train_size=None, test_size=None, random_state=None, shuffle=True):
    '''returns two kernel lists, for train and test'''

    idx = range(len(Y))
    train,test = skms.train_test_split(idx,
        train_size=train_size,
        test_size=test_size,
        random_state=random_state,
        shuffle=shuffle)
    KL_tr = [K[train][:,train] for K in KL]
    KL_te = [K[test ][:,train] for K in KL]
    return KL_tr, KL_te, Y[train], Y[test]
