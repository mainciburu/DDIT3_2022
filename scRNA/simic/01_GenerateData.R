##########################################
#######   Prepare data for Simic   #######
##########################################

library(reticulate)
reticulate::use_python("/home/mainciburu/scRNA/simic/myenv/bin/python3")
py_discover_config("magic")
library(Rmagic)

library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)


seurat.ery <- readRDS('/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds')

# Analysis of the data: DDIT3 vs control
cell_populations <- c('DDIT3', 'control')

# Remove reticulocytes from the analysis
Idents(seurat.ery) <- seurat.ery$Celltype
seurat.ery <- subset(seurat.ery, idents=c("CD34", "BFU", "CFU", "Proerythroblast", "Early_basophilic", "Late_basophilic", "Polychromatic", "Orthochromatic"))
Idents(seurat.ery) <- seurat.ery$Condition

# keep only the labels with more than 10 representants in both groups
# table(seurat.ery$labels, seurat.ery$orig.ident)
stats <- setNames(as.data.frame(table(seurat.ery$singleR[seurat.ery$Condition == 'control']), stringsAsFactors =F), c('singleR', 'control'))
stats <- merge( stats, setNames(as.data.frame(table(seurat.ery$singleR[seurat.ery$Condition == 'DDIT3']), stringsAsFactors =F), c('singleR', 'DDIT3')), by= 'singleR')
labels_2_keep <- stats[stats$control > 10 & stats$DDIT3 >10,'singleR']

cells_2_keep <- names(seurat.ery$singleR[seurat.ery$singleR %in% labels_2_keep])

### data from seurat ####
seurat.ery_raw <- as.data.frame(seurat.ery@assays$RNA@counts)
# remove unexpressed genes
seurat.ery_raw <- seurat.ery_raw[, colnames(seurat.ery_raw) %in% cells_2_keep]
unexpresed_genes <- names(which(rowSums(abs(seurat.ery_raw))<1e-6))
seurat.ery_raw <- seurat.ery_raw[ !rownames(seurat.ery_raw) %in% unexpresed_genes, ] #cells_2_keep]
# remove undesired cells

# cell assignments
cell_group_idx <-  as.data.frame(seurat.ery$Condition)
cell_group_idx <- as.character(cell_group_idx[rownames(cell_group_idx) %in% cells_2_keep,])
assignment <- as.character(seq(0,length(cell_populations)-1))
names(assignment) <- cell_populations
cluster_assignments <- assignment[cell_group_idx]


### RUN MAGIC 
seurat.ery_raw<-Rmagic::library.size.normalize(t(seurat.ery_raw))
seurat.ery_raw <- sqrt(seurat.ery_raw)
data_MAGIC <- magic(seurat.ery_raw,genes='all_genes') 

data_MAGIC_df <- as.data.frame(data_MAGIC)

data_MAGIC_df <- data_MAGIC_df

analysis_root_dir <- '/home/mainciburu/scRNA/simic/data/'
file_idx <- 'simic_ddit3_filtered'
DF_f <- paste0(analysis_root_dir,file_idx,".DF.pickle") #a genes x cells matrix
TFs_f <- paste0(analysis_root_dir,file_idx,".TFs.pickle") #name of genes to be use as drivers
cluster_assignment_f <- paste0(analysis_root_dir,file_idx,".clustAssign.txt") #phenotype of interes of the cells


### SELECT TFs and TARGETS ####
MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000

TFs_list <- py_load_object('/home/mainciburu/scRNA/simic/data//human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

# select top100 most variable TFs
MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)    
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

# select top1000 most variable targets (exclude TFs)
target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

sum(is.na(TFs))
### WRITE DATA TO DISK ####

write(cluster_assignments,file = cluster_assignment_f)

sum(is.nan(as.matrix(input_data)))
sum(is.infinite(as.matrix(input_data)))

reticulate::py_save_object(as.data.frame(input_data), filename = DF_f)
reticulate::py_save_object(TFs,filename = TFs_f)


