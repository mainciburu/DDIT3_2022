#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 jianhao2 <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.

############
# 1) Run simicLASSO with cross_val=False
# 2) Choose lambdas => lambda1 = 0.01, lambda2 = 0.1
# 3) Run simicLASSO with lambda1 = 0.01, lambda2 = 0.1 cross_val=False
# 4) Filter weights in R (Filter_weigths.R)
# 5) Calculate AUC from filtered weigths (ddit3_weights_filtered_BIC.pickle)
# 6) AUC pickle to txt 
# 7) Analyze AUC in R
#############

from simiclasso.clus_regression import simicLASSO_op
from simiclasso.weighted_AUC_mat import main_fn
import time

p2df = '/home/mainciburu/scRNA/simic/data/simic_ddit3_filtered.DF.pickle'
p2assignment = '/home/mainciburu/scRNA/simic/data/simic_ddit3_filtered.clustAssign.txt'
k_cluster = 2
similarity = True
p2tf = '/home/mainciburu/scRNA/simic/data/human_TF.pickle'
p2saved_file = '/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_weights.pickle'
num_TFs = 100
num_target_genes = 1000
max_rcd_iter = 50000
df_with_label = False
_NF = 100          # normalization factor

# simicLASSO_op(p2df, p2assignment, k_cluster, similarity, p2tf, p2saved_file, 
#              num_TFs, num_target_genes, _NF = _NF, 
#              max_rcd_iter = max_rcd_iter, df_with_label = df_with_label)
# 


lambda1=0.01
lambda2=0.1

simicLASSO_op(p2df, p2assignment, similarity, p2tf, p2saved_file,  
			k_cluster, num_TFs, num_target_genes, 
			max_rcd_iter = max_rcd_iter, df_with_label = False, 
			lambda1=lambda1, lambda2 = lambda2, cross_val=False)



p2df = '/home/mainciburu/scRNA/simic/data/simic_ddit3_filtered.DF.pickle'
p2AUC = '/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_AUCs_filtered_BIC.pickle'
p2saved_file = "/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_weights_filtered_BIC.pickle"
main_fn(p2df, p2saved_file, p2AUC)


