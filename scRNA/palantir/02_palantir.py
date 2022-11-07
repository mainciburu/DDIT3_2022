import sys
sys.path.insert(0, "/home/mainciburu/scRNA/scanpy/venv/lib/python3.6/site-packages/")

import numpy as np
import pandas as pd
import pickle
import matplotlib
import matplotlib.pyplot as plt
import scipy
import statsmodels
import palantir
import scanpy as sc
import seaborn as sns
from collections import Counter
import feather

########## Input data ####################

### DDIT3 dataset
h5ad_file="/home/mainciburu/scRNA/diff_ddit3/loom/ddit3_ery.h5ad"
loom_file="/home/mainciburu/scRNA/diff_ddit3/loom/ddit3_ery_norm.loom"
norm_df_file="/home/mainciburu/scRNA/diff_ddit3/loom/ddit3_ery_norm_mat.csv"
norm_df_ori_file="/home/mainciburu/scRNA/diff_ddit3/loom/ddit3_ery_norm_mat_full.csv"
umap_seurat_file="/home/mainciburu/scRNA/diff_ddit3/umap_ddit3_ery.txt"

plot_name = "/home/mainciburu/scRNA/diff_ddit3/pics/palantir/ddit3_ery"
res_path="/home/mainciburu/scRNA/diff_ddit3/palantir_results/ddit3/"

### Control dataset
#h5ad_file="/home/mainciburu/scRNA/diff_ddit3/loom/control_ery.h5ad"
#loom_file="/home/mainciburu/scRNA/diff_ddit3/loom/control_ery_norm.loom"
#norm_df_file="/home/mainciburu/scRNA/diff_ddit3/loom/control_ery_norm_mat.csv"
#norm_df_ori_file="/home/mainciburu/scRNA/diff_ddit3/loom/control_ery_norm_mat_full.csv"
#umap_seurat_file="/home/mainciburu/scRNA/diff_ddit3/umap_control_ery.txt"

#plot_name = "/home/mainciburu/scRNA/diff_ddit3/pics/palantir/control_ery"
#res_path="/home/mainciburu/scRNA/diff_ddit3/palantir_results/control/"

##########################################

# load scanpy object processed
adata = sc.read_loom(loom_file)
adata.write(filename=h5ad_file)
adata = sc.read(filename=h5ad_file)

# Integrated matrix from csv
norm_df=pd.read_csv(norm_df_file, index_col=0)
# transpose
norm_df=norm_df.T

# Project results on seurat umap
umap_seurat=pd.read_csv(umap_seurat_file, sep="\t", header=None, index_col = 0)
umap_seurat.columns = ['x', 'y']

# Preprocess data for Palantir
# PCA
pca_projections, _ = palantir.utils.run_pca(norm_df)
# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=20)
# Multiscale data
ms_data = palantir.utils.determine_multiscale_space(dm_res)

start_cell=["TTGCTGCGTTCTCAGA_DDIT3"]         ### DDIT3
#start_cell=["ATCGCCTGTCGTGATT_control"]      ### Control

pr_res = palantir.core.run_palantir(ms_data, start_cell, n_jobs = 1)

final_cell=list(pr_res.branch_probs.columns)

# Rename trajectories
i=pr_res.branch_probs.columns
df = adata.obs['singleR']
df = df.loc[i]
names = list(df)
pr_res.branch_probs.columns = names

# Plot results
cells=[start_cell, final_cell]
palantir.plot.highlight_cells_on_tsne(umap_seurat, final_cell)
plt.savefig(plot_name + "_umap_final_cells.pdf")

palantir.plot.plot_palantir_results(pr_res, umap_seurat)
plt.savefig(plot_name + "_umap_branches.pdf")

# pseudotime, branch probability, differentiation potential
pr_res.branch_probs.to_csv(res_path + "branch_probs.csv")
pr_res.pseudotime.to_csv(res_path + "pseudotime.csv")
pr_res.entropy.to_csv(res_path + "diff_potential.csv")

# Save pr_res object
file_pr = open(res_path + 'pr_res.obj', 'wb') 
pickle.dump(pr_res, file_pr)

# Gene expression trends
# !!Use original normalized expression values (not integrated)
norm_df_ori = pd.read_csv(norm_df_ori_file, index_col = 0)
norm_df_ori = pd.DataFrame.transpose(norm_df_ori)

# Impute data
imp_df = palantir.utils.run_magic_imputation(norm_df_ori, dm_res)
feather.write_dataframe(imp_df, res_path + 'imp_df.feather')

