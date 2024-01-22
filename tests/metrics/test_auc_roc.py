import sys

sys.path.append('/projectnb2/docking/imhaoyu/AbEMap-explore')
import csv
from tools import seq_align, metrics
import glob

sys.path.append('/projectnb2/docking/imhaoyu/my_protein_toolbox')
from metrics import classification
import numpy as np
import unittest

native_pdb = '/projectnb2/docking/imhaoyu/AbEMap-explore/draft_AbEMap_data/BM5_antigen_native/1BVK_l_b.pdb'
abemap_pdb = '/projectnb2/docking/imhaoyu/AbEMap-explore/draft_AbEMap_data/BM5_AbEMap/1BVK.pdb'

AbE_structure = seq_align.get_structure_pdb(abemap_pdb)
native_structure = seq_align.get_structure_pdb(native_pdb)
bfactor_AbE, bfactor_true = seq_align.align_seq_bfactor(AbE_structure, native_structure, sum=False)
y_score = []
y_true = []
for t1, t2 in zip(bfactor_AbE, bfactor_true):
    y_score.append(t1[1])
    y_true.append(t2[1])

classification._plot_roc_auc(np.array(y_true), np.array(y_score))
print(classification._roc_auc_score(np.array(y_true), np.array(y_score)))

print(classification.sklearn_roc_auc(y_true, y_score))
# print(metrics.cal_confusion_matrix_scores(bfactor_AbE, bfactor_true, top=None))
