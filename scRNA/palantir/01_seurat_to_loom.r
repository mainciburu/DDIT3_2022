# Extract from Seurat object
  # loom object
  # Normalized expression matrix
  # umap coordinates

library(loomR)
library(Seurat)
#library(hdf5r)

###
seurat.file<-"/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds"
seurat<-readRDS(seurat.file)

res.path<-"/home/mainciburu/scRNA/diff_ddit3/loom/"
res.name<-"ddit3"
#res.name<-"control"

####### .loom objects #########
# Remove meta.data columns with NA
seurat@meta.data<-seurat@meta.data[,colSums(is.na(seurat@meta.data))==0]
Idents(seurat)<-"Condition"
seurat<-subset(seurat, idents = "DDIT3")
#seurat<-subset(seurat, idents = "Control")

# Create .loom object with original normalized expression matrix 
dat<-CreateSeuratObject(counts = seurat@assays$RNA@data, assay = "RNA", meta.data = seurat@meta.data)
dat@assays$RNA@data<-dat@assays$RNA@counts
dat@assays$RNA@var.features<-rownames(dat@assays$RNA@counts)
dat.loom<-as.loom(dat, filename = paste0(res.path, res.name, "_norm.loom"))
dat.loom$close_all()

####### umap ###########
# Get integrated umap coordinates
umap<-seurat@reductions$umap.int@cell.embeddings
write.table(x = umap, file = paste0(res.path, res.name, "_umap.txt"), sep = "\t", quote = F, col.names = F, row.names = T)


######## csv expression matrix #############

# Variable genes normalized matrix
x<-seurat@assays$RNA@data[seurat@assays$RNA@var.features,]
write.csv(x = x, file = paste0(res.path, res.name, "_norm_mat.csv"))

# Full original normalized matrix
x<-seurat@assays$RNA@data
write.csv(x = x, file = paste0(res.path, res.name, "_norm_mat_full.csv"))