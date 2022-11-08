########################################
## Preprocess individual samples
## Integrate
## SingleR predictions
## Reintegrate erythroid data
########################################

library(Seurat)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(umap)
library(dplyr)
library(cluster)
library(future)

col.ery<-c("CD34"="#FFFF33", "BFU"="#A65628",
           "CFU"="#F781BF", "Proerythroblast"="#E41A1C", 
           "Early_basophilic"="#377EB8", "Late_basophilic"="#4DAF4A", 
           "Polychromatic"="#FF7F00", "Orthochromatic"="#984EA3", "Reticulocytes"="#899DA4")

############ FUNCTIONS ####################################

# init_seurat - create seurat object and calculate % of mitochondrial genes
# input
	# name: sample name
	# path: path to the filtered 10X matrix
	# cols: string to add to the barcode names
	# condition: Young, Senior or MDS
# output: seurat object 

init_seurat<-function(name, path, cols, condition){
	init.data<-Read10X(data.dir = path)
	colnames(init.data)<-paste0(colnames(init.data), "_", cols)
	seurat.obj<-CreateSeuratObject(counts = init.data, min.cells = 3, min.features = 100, project = name)
	seurat.obj$Patient<-name
	seurat.obj$Condition<-condition

	mito.genes<-grep("MT-", rownames(seurat.obj), value = T)
	percent.mito<-Matrix::colSums(GetAssayData(seurat.obj, "counts")[mito.genes,])/Matrix::colSums(seurat.obj)
	seurat.obj$percent.mito<-percent.mito
	return(seurat.obj)
}

# preprocess_seurat - normalize, calculate cell cycle scores, scale, regress effects, find variable genes and PCA
# input
	# seurat.obj: object
	# npcs: number of principal component to calculate
	# vars: variables from metadata to regress
# output: processed seurat object

preprocess_seurat<-function(seurat.obj, npcs, vars){
	seurat.obj<-NormalizeData(seurat.obj)
	seurat.obj<-CellCycleScoring(object = seurat.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
	seurat.obj<-ScaleData(object = seurat.obj, vars.to.regress = vars)
	seurat.obj<-FindVariableFeatures(seurat.obj, selection.method = "vst",
                            nfeatures = 2000)
	seurat.obj<-RunPCA(seurat.obj, npcs = npcs)
	return(seurat.obj)
}

# cluster_seurat: cluster cells with specific resolution, plot and calculate average silhouette
# input
	# seurat.obj: object
	# npcs: number of significant components
	# resolution: clustering resolution
# output: seurat object with clusters, cluster umap plot and silhouette barplot

cluster_seurat<-function(seurat.obj, npcs, resolution){
	seurat.obj<-FindNeighbors(object = seurat.obj, reduction = "pca", dims = 1:npcs)
	seurat.obj<-FindClusters(seurat.obj, resolution = resolution)
	col<-colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat.obj$seurat_clusters)))
	print(DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", cols = col))
	# Cluster silhouette
	d<-dist(x = seurat.obj@reductions$pca@cell.embeddings[,1:npcs])
	s<-silhouette(x = as.numeric(seurat.obj$seurat_clusters), dist = d)
	summary(s)
	s.avg<-as.numeric(summary(s)$clus.avg.widths)
	c<-length(unique(seurat.obj$seurat_clusters)) - 1
	barplot(s.avg, horiz = T, names.arg = as.character(0:c), col = col)
	return(seurat.obj)
}


####################################################

###############################
#######   Preprocess  #########
####### DDIT3 dataset #########
###############################

## INPUT 
name<-"DDIT3"
path<-"/home/mainciburu/data/SC_MDS/DDIT3/diff_DDIT3_Count/outs/filtered_feature_bc_matrix/"
cols<-"DDIT3"
condition<-"DDIT3"

# Create seurat object
seurat<-init_seurat(name, path, cols, condition)

# Filter 
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mito", feature2 = "nFeature_RNA")
VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))

mito.low<-0
mito.high<-0.15
count.low<-0
count.high<-60000    
feature.low<-0
feature.high<-Inf
seurat<-subset(seurat, percent.mito>mito.low & percent.mito<mito.high & 
			   nCount_RNA > count.low & nCount_RNA < count.high &
			   nFeature_RNA > feature.low & nFeature_RNA < feature.high)

# Activate parallelization
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 2097152000)

# Normalize, scale, regress effects and PCA
npcs<-30
vars<-c("S.Score", "G2M.Score")
seurat<-preprocess_seurat(seurat, npcs, vars)

# Decide dimmensionality
ElbowPlot(seurat, ndims = 30)
DimPlot(seurat, reduction = "pca", group.by = "Phase")

# umap
npcs<-15
seurat<-RunUMAP(object = seurat, assay = "RNA", dims = 1:npcs)

# explore unwanted effects
DimPlot(seurat, reduction = "umap", group.by = "Phase")
FeaturePlot(seurat, 
            features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
            reduction = "umap", cols = rev(brewer.pal(11, "RdYlBu")))

# Save object
saveRDS(seurat, file = "seurat_ddit3.rds")


###############################
#######   Preprocess  #########
####### control dataset #######
###############################

## INPUT 
name<-"control"
path<-"/home/mainciburu/data/SC_MDS/DDIT3/diff_control_Count/outs/filtered_feature_bc_matrix/"
cols<-"control"
condition<-"control"

# Create seurat object
seurat<-init_seurat(name, path, cols, condition)

# Filter 
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat, feature1 = "percent.mito", feature2 = "nFeature_RNA")
VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))

mito.low<-0
mito.high<-0.15
count.low<-0
count.high<-75000  
feature.low<-0
feature.high<-Inf
seurat<-subset(seurat, percent.mito>mito.low & percent.mito<mito.high & 
			   nCount_RNA > count.low & nCount_RNA < count.high &
			   nFeature_RNA > feature.low & nFeature_RNA < feature.high)

# Activate parallelization
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 2097152000)

# Normalize, scale, regress effects and PCA
npcs<-30
vars<-c("S.Score", "G2M.Score")
seurat<-preprocess_seurat(seurat, npcs, vars)

# Decide dimmensionality
ElbowPlot(seurat, ndims = 30)
DimPlot(seurat, reduction = "pca", group.by = "Phase")

# umap
npcs<-15
seurat<-RunUMAP(object = seurat, assay = "RNA", dims = 1:npcs)

# explore unwanted effects
DimPlot(seurat, reduction = "umap", group.by = "Phase")
FeaturePlot(seurat, 
            features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
            reduction = "umap", cols = rev(brewer.pal(11, "RdYlBu")))

# Save object
saveRDS(seurat, file = "seurat_control.rds")



###############################
#######  Integration  #########
###############################

seurat1<-readRDS("seurat_control.rds")
seurat2<-readRDS("seurat_ddit3.rds")

# List of objects
AllData<-list(seurat1, seurat2)
names(AllData)<-c("seurat1", "seurat2")

# Activate parallelization
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 2097152000)     ## change

# integration
nfeatures<-2000
ndim<-50
int.features<-SelectIntegrationFeatures(object.list = AllData, nfeatures = nfeatures)

anchors<-FindIntegrationAnchors(object.list = AllData, 
                                dims = 1:ndim)
seurat.int<-IntegrateData(anchorset = anchors, dims = 1:ndim)


# Analysis on integrated data

# Recalculate cell Cycle score 
DefaultAssay(seurat.int)<-"RNA"
seurat.int<-CellCycleScoring(seurat.int, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# scale integrated data and regress cc effect
DefaultAssay(seurat.int)<-"integrated"
seurat.int<-ScaleData(seurat.int, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "nFeature_RNA"))

# PCA
seurat.int<-RunPCA(seurat.int, npcs = ndim)
DimPlot(seurat.int, dims = 1:2, reduction = "pca", group.by ="Condition", pt.size = 0.1)
# decide dimensionality
ElbowPlot(seurat.int, ndims = ndim)

# umap
npcs<-20
reduction.name<-"umap"
seurat.int<-RunUMAP(object = seurat.int, assay = "integrated", dims = 1:npcs, reduction.name = reduction.name)
DimPlot(seurat.int, reduction = "umap", dims = 1:2, group.by = "Condition", pt.size = 0.1)

# clustering  
resolution<-0.4
seurat.int<-cluster_seurat(seurat.int, npcs, resolution)

# markers
cl.markers<-FindAllMarkers(seurat.int, assay = "RNA", test.use = "MAST", only.pos = T)
top.markers<-cl.markers%>%group_by(cluster)%>%top_n(avg_logFC, n = 5)
condition.markers<-FindMarkers(object = seurat.int, ident.1 = "DDIT3", ident.2 = "control", 
                               test.use = "MAST", only.pos = F)


# Ribosomal genes
rps<-grep(pattern = "^RPS", x = rownames(seurat.int))
rpl<-grep(pattern = "^RPL", x = rownames(seurat.int))
ribo.genes<-c(rps, rpl)
percent.ribo<-Matrix::colSums(GetAssayData(seurat.int, "counts")[ribo.genes,])/Matrix::colSums(GetAssayData(seurat.int, "counts"))
seurat.int$percent.ribo<-percent.ribo

# save object
saveRDS(seurat.int, file = "seurat_int.rds")


##########################################
#######  SingleR Classification  #########
##########################################

library(SingleR)
library(DESeq2)

# Reference
load("scRNA/diff_ddit3/counts_rnaseq_pb.Rdata")
dds<-DESeqDataSetFromMatrix(countData = counts, colData = info, design = ~state)
dds<-dds[rowSums(assay(dds))!=0]
dds<-DESeq(dds, test = "LRT", reduced = ~1)
ref<-rlog(dds)
ref<-assay(ref)                      

labs<-dds$state

# Test
seurat<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int.rds")
test<-seurat@assays$RNA@data
test<-as.matrix(test)

# singleR
pred<-SingleR(test = test, ref = ref, labels = labs)

seurat$singleR<-pred$labels
seurat$singleR<-factor(seurat$singleR, levels = unique(dds$state))
seurat$singleR.score<-pred$tuning.scores$first

# heatmap scores
annot<-list(Labels=col.ery)
plotScoreHeatmap(pred, treeheight_row = 0, annotation_colors = annot)

p1<-DimPlot(seurat, group.by = "singleR", reduction = "umap", cols = col.ery)
p2<-FeaturePlot(seurat, features = "singleR.score", cols = rev(brewer.pal(11, "RdYlBu")))
plot_grid(p1, p2, nrow = 1)

# Markers
Idents(seurat)<-"singleR"
cl.markers<-FindAllMarkers(seurat, only.pos = T, assay = "RNA")
top.10<-cl.markers%>%group_by(cluster)%>%top_n(avg_logFC, n = 10)

# condition deg
seurat$singleR_condition<-paste0(as.character(seurat$singleR), 
                                 "_", seurat$Condition)
Idents(seurat)<-"singleR_condition"
DefaultAssay(seurat)<-"RNA"
res<-data.frame()
for(state in unique(seurat$singleR)){
  res.i<-FindMarkers(seurat, ident.1 = paste0(state, "_DDIT3"), 
                   ident.2 = paste0(state, "_control"), test.use = "MAST", 
                   logfc.threshold = 0, min.pct = 0, assay = "RNA")
  #res.i<-res.i[res.i$p_val_adj<0.05,]
  res.i$cluster<-state
  res.i$gene<-rownames(res.i)
  colnames(res.i)[c(3,4)]<-c("pct.DDIT3", "pct.control")
  res<-rbind(res, res.i)
}

saveRDS(seurat, "/home/mainciburu/scRNA/diff_ddit3/seurat_int.rds")


##########################################
#######  Erythroid clusters only #########
##########################################


## Delete myeloid clusters: 7 and 9 
# 7 => mix of myeloid progenitors
# 9 => basophil progenitors
Idents(int)<-"integrated_snn_res.0.4"
int.ery<-subset(int, idents = c(0:6,8))

# Reintegrate by condition
AllData<-SplitObject(object = int.ery, split.by = "Condition")

# Normalize, scale, find variable genes
for(i in 1:2){
  AllData[[i]]<-NormalizeData(AllData[[i]], scale.factor = 10000)
  AllData[[i]]<-ScaleData(AllData[[i]], features = rownames(AllData[[i]]))
  AllData[[i]]<-FindVariableFeatures(AllData[[i]], selection.method = "vst", nfeatures = 3000)
}

# Integrate
anchors<-FindIntegrationAnchors(AllData, dims = 1:30)
int.ery<-IntegrateData(anchors, features.to.integrate = rownames(AllData$control), dims = 1:30)

# Scale integrated data and regress cell cycle and count effects
int.ery@active.assay<-"RNA"
int.ery<-CellCycleScoring(int.ery, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
int.ery@active.assay<-"integrated"
int.ery<-ScaleData(object = int.ery, vars.to.regress = c("S.Score", "G2M.Score"))

# UMAP
int.ery<-RunPCA(int.ery, npcs = 30)
ElbowPlot(object = int.ery, ndims = 30)
int.ery<-RunUMAP(int.ery, dims = 1:15, seed.use = 1234)
DimPlot(int.ery, reduction = "umap", group.by = "Condition")
DimPlot(int.ery, reduction = "umap", group.by = "Phase")
DimPlot(int.ery, reduction = "umap", group.by = "singleR")
DimPlot(int.ery, reduction = "umap", group.by = "integrated_snn_res.0.4")

FeaturePlot(int.ery, reduction = "umap", features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"))

## Combine singleR and reticulocytes cluster
int.ery$CellType<-int.ery$singleR
int.ery$CellType[int.ery$integrated_snn_res.0.4==4]<-"Reticulocytes"
int.ery$CellType<-factor(int.ery$CellType, levels = c("CD34", "BFU", "CFU", "Proerythroblast", "Early_basophilic", "Late_basophilic", "Polychromatic", "Orthochromatic", "Reticulocytes"))

saveRDS(int.ery, file = "scRNA/diff_ddit3/seurat_int_ery_reti.rds")

# differential expression per condition and cluster
DefaultAssay(int.ery)<-"RNA"
Idents(int.ery)<-"singleR_condition"
res<-data.frame()
for(state in unique(int.ery$singleR)){
  res.i<-FindMarkers(int.ery, ident.1 = paste0(state, "_DDIT3"),
                     ident.2 = paste0(state, "_control"), test.use = "MAST",
                     logfc.threshold = 0, min.pct = 0)
  res.i$cluster<-state
  res.i$gene<-rownames(res.i)
  colnames(res.i)[c(3,4)]<-c("pct.DDIT3", "pct.control")
  res<-rbind(res, res.i)
}
saveRDS(res, file = "scRNA/diff_ddit3/deg_condition_ery_clusters_all.rds")








