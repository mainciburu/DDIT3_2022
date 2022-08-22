##################################
#### Generate data for scVelo ####
##################################


library(loomR)
library(Seurat)

seurat<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds")

# Remove meta.data columns with NA
seurat@meta.data<-seurat@meta.data[,colSums(is.na(seurat@meta.data))==0]

### Loom object
# Both conditions => separate later
## RNA matrix
dat<-CreateSeuratObject(counts = seurat@assays$RNA@data, assay = "RNA", meta.data = seurat@meta.data)
dat@assays$RNA@data<-dat@assays$RNA@counts
dat@assays$RNA@var.features<-rownames(dat@assays$RNA@counts)
dat.loom<-as.loom(dat, filename = "/home/mainciburu/scRNA/diff_ddit3/loom/ery_reti_RNA.loom")
dat.loom$close_all()

## Seurat UMAP
umap<-seurat@reductions$umap@cell.embeddings
rownames(umap)<-gsub("_control", "", rownames(umap))
rownames(umap)<-gsub("_DDIT3", "", rownames(umap))
write.table(x = umap, file = "/home/mainciburu/scRNA/diff_ddit3/scvelo/ery_reti_umap_scv.txt", sep = "\t", quote = F, col.names = F, row.names = T)
