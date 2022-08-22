######## scvelo ############

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy
import scanpy as sc
import scvelo as scv
import loompy

########## Input data ####################
adata_file="/home/mainciburu/scRNA/diff_ddit3/loom/ery_reti_RNA.loom"    # seurat processed object
scvdata_file="/home/mainciburu/scRNA/diff_ddit3/loom/control_ddit3_Count.loom"    # velocyto control and ddit3 combined matrices
umap_seurat_file="/home/mainciburu/scRNA/diff_ddit3/ery_reti_umap.txt"      # umap coordinates

plot_name = "/home/mainciburu/scRNA/diff_ddit3/pics_velocyto/3.0/"
res_path="/home/mainciburu/scRNA/diff_ddit3/scvelo/"

col_control=["#A65628", "#F781BF", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C", "#899DA4"]
col_ddit3=["#A65628", "#FFFF33", "#F781BF", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C", "#899DA4"]
###########################################

######### store and read results ##########
#adata.write(res_path + 'ery_reti_scv.h5ad', compression='gzip')
#adata = scv.read(res_path + 'ery_reti_scv.h5ad')
#adata_ddit3 = scv.read(res_path + 'ddit3_ery_reti_scv.h5ad')
#adata_control = scv.read(res_path + 'control_ery_reti_scv.h5ad')

# multiple scvdata => merge them
files = ["/home/mainciburu/data/SC_MDS/DDIT3/diff_control_Count/velocyto/diff_control_Count.loom",
		 "/home/mainciburu/data/SC_MDS/DDIT3/diff_DDIT3_Count/velocyto/diff_DDIT3_Count.loom"]
loompy.combine(files, scvdata_file, key="Accession")

####################################
######### Preprocessing ############
####################################

# Settings
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# Read data
adata = sc.read(adata_file)
scvdata = scv.read(scvdata_file, cache=True)
adata = scv.utils.merge(adata, scvdata)

# umap from seurat
umap_seurat=pd.read_csv(umap_seurat_file, sep="\t", header=None, index_col = 0)
umap_seurat.index = adata.obs.index
#umap_seurat = umap_seurat.loc[adata.obs.index,]
umap_seurat = umap_seurat.to_numpy()
adata.obsm["X_umap"] = umap_seurat

# Plot proportion of unspliced RNA
counts_s = scv.utils.sum_var(adata.layers['spliced'])
counts_u = scv.utils.sum_var(adata.layers['unspliced'])
fractions_u = counts_u / (counts_s + counts_u)
scv.pl.scatter(adata, basis="umap", color=fractions_u, figsize = (10,5), size=80)
plt.savefig(plot_name + "prop_unspliced.pdf")
scv.pl.proportions(adata, groupby="singleR", figsize = (20,5))                                                                                                                                          
plt.savefig(plot_name + "_prop_unspliced_summary.pdf")  

# Preprocessing
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)  # ***already normalized in seurat
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Split in two objects
adata_control = adata[adata.obs['Condition'] == 'control']
adata_ddit3 = adata[adata.obs['Condition'] == 'DDIT3']

# Velocity dynamical model - independently per condition
scv.pp.neighbors(adata_control, n_neighbors = 30)
scv.tl.recover_dynamics(adata_control, var_names="velocity_genes")      # infer full kinetics (transcription, degradation, splicing rates). Needed for dynamical model
										                                # fits and velocity_genes are established here
scv.tl.velocity(adata_control, mode='dynamical')
scv.tl.velocity_graph(adata_control)    

scv.pp.neighbors(adata_ddit3, n_neighbors = 30)
scv.tl.recover_dynamics(adata_ddit3, var_names="velocity_genes")      
scv.tl.velocity(adata_ddit3, mode='dynamical')
scv.tl.velocity_graph(adata_ddit3)    

adata_control.write(res_path + 'control_ery_reti_scv.h5ad', compression='gzip')
adata_ddit3.write(res_path + 'ddit3_ery_reti_scv.h5ad', compression='gzip')


###############################################################################

# project velocities on UMAP
scv.pl.velocity_embedding_stream(adata_control, basis='umap', color="CellType", alpha = 0.5, 
	legend_loc='upper left', fontsize = 12, figsize = (8,5), min_mass=1, palette=col_control,
	save = "control_velocity_stream.png", dpi=1000)
scv.pl.velocity_embedding_stream(adata_ddit3, basis='umap', color="CellType", alpha = 0.5, 
	legend_loc='upper left', fontsize = 12, figsize = (8,5), min_mass=1, palette=col_ddit3,
	save = "ddit3_velocity_stream.png", dpi=1000)


scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, basis='umap', color="singleR")
plt.savefig(plot_name + "_velocity_scatter.pdf")


# length of velocity -> speed or rate of differentiation
# confidence -> coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities)
scv.tl.velocity_confidence(adata_control)
scv.tl.velocity_confidence(adata_ddit3)

keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata_control, c=keys, cmap='coolwarm', perc=[5, 95])
scv.pl.scatter(adata_ddit3, c=keys, cmap='coolwarm', perc=[5, 95])

plt.savefig(plot_name + "_speed_coherence.pdf")
df = adata.obs.groupby('singleR')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

# store results 
adata.write(res_path + 'loom/control_ery_scv.h5ad', compression='gzip')

