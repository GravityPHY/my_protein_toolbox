import numpy as np
import matplotlib.pyplot as plt
from statistics import mode


def F1(TP, FP, FN):
    return TP / (TP + 0.5 * (FP + FN))


def MCC(TP, TN, FP, FN):
    try:
        _mcc=(TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        return _mcc
    except:
        return 0


def find_top_cutoff(tuple, top=10):
    """
    find the cutoff for the first top values in the list,
    if top was set to a value longer than the length of tuple, it will be the minimum in the list
    (in order to avoid out of index error)
    :param tuple:
    :param top:
    :return:
    """
    tuple_length = len(tuple)
    new_tuple = sorted(tuple, key=lambda x: x[1], reverse=True)
    if top < tuple_length:
        return new_tuple[top][1]
    else:
        return new_tuple[-1][1]


def cal_confusion_matrix_scores(bfactor_pred, bfactor_true, top=None):
    try:
        assert len(bfactor_pred) == len(bfactor_true)
    except:
        print(len(bfactor_pred), '!=', len(bfactor_true))
    TP, TN, FP, FN = 0, 0, 0, 0
    tot = 0
    cutoff = 0.1
    if top is not None:
        cutoff = find_top_cutoff(bfactor_pred, top)
    for t1, t2 in zip(bfactor_pred, bfactor_true):
        if t2[0] != '0':
            # tot += 1
            if t1[1] >= cutoff and int(t2[1]) == 1:
                TP += 1
            elif t1[1] >= cutoff and int(t2[1]) == 0:
                FP += 1
            elif int(t2[1]) == 0:
                TN += 1
            elif int(t2[1]) == 1:
                FN += 1
    # assert tot == TP + TN + FP + FN
    return TP, FP, FN, TN


def top_confusion_matrix(bfactor_1, bfactor_2, bfactor_true, top1, top2):
    """
    calculate the confusion matrix for top1 residues in structure1 and
    top2 residues in structure 2
    :param bfactor_1: (seq, bfactor) list from structure1
    :param bfactor_2: (seq, bfactor) list form structure2
    :param bfactor_true: (seq, bfactor) list from true structure
    :int top1: number of top ranking residues from structure 1, must be an integer
    :int top2: number of top ranking residues from structure 2, must be an integer
    :return: TP - true positive, TN - true negative, FP - false positive, FN - false negative
    """
    assert isinstance(top1, int)
    assert isinstance(top2, int)
    TP, TN, FP, FN = 0, 0, 0, 0
    tot = 0
    cutoff1 = find_top_cutoff(bfactor_1, top1)
    cutoff2 = find_top_cutoff(bfactor_2, top2)
    for t1, t2, t in zip(bfactor_1, bfactor_2, bfactor_true):
        if t[0] != '0':
            #assert t1[0] == t2[0]
            if t1[1] >= cutoff1 and t2[1] >= cutoff2:
                if int(t[1]) == 1:
                    TP += 1
                if int(t[1]) == 0:
                    FP += 1
            else:
                if int(t[1]) == 1:
                    FN += 1
                elif int(t[1]) == 0:
                    TN += 1
            # if t1[1] >= cutoff1 and t2[1] >= cutoff2 and int(t[1]) == 1:
            #    TP += 1
            # elif t1[1] >= cutoff1 and t2[1] >=cutoff2 and int(t[1]) == 0:
            #    FP += 1
            # elif t1[1] >= cutoff1 and t2[1] < cutoff2 and int(t[1]) == 1:
            #    FN += 1
            # elif t1[1] < cutoff1 and t2[1] >= cutoff2 and int(t[1]) == 1:
            #    FN += 1
            # elif t1[1] >= cutoff1 and t2[1] < cutoff2 and int(t[1]) == 0:
            #    TN += 1
            # elif t1[1] < cutoff1 and t2[1] >= cutoff2 and int(t[1]) == 0:
            #    TN += 1
    return TP, FP, FN, TN
