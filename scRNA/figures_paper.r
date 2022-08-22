######## Figures ##########
library(Seurat)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(xlsx)
library(plyr)
col.ery<-c("CD34"="#FFFF33", "BFU"="#A65628",
           "CFU"="#F781BF", "Proerythroblast"="#E41A1C", 
           "Early_basophilic"="#377EB8", "Late_basophilic"="#4DAF4A", 
           "Polychromatic"="#FF7F00", "Orthochromatic"="#984EA3", "Reticulocytes"="#899DA4")

source("/home/mainciburu/scRNA/pipeline_2.0/08_trajectory_analysis/plot_gene_trends.r")


########################
#####   Figure 4   #####
########################

# UMAP with labels includint reticulocytes
int.ery.reti<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds")
pdf("scRNA/diff_ddit3/figures_paper/umap_labels.pdf", 
    useDingbats = F, width = 15, height = 8)
pp<-DimPlot(int.ery.reti, group.by = "CellType", cols = col.ery, pt.size = 0.8)
pp<-pp + theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
  theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
print(pp)
dev.off()

# Barplot with proportion per identity 
tt<-prop.table(table(int.ery.reti$CellType, int.ery.reti$Condition), margin = 2)
df<-melt(tt)
colnames(df)<-c("Identity", "Condition", "Proportion")
df$Identity<-factor(df$Identity, levels = levels(int.ery.reti$CellType))
df$Condition<-factor(df$Condition, levels = c("control", "DDIT3"))
pdf("scRNA/diff_ddit3/figures_paper/proportions.pdf", 
    useDingbats = F, width = 12, height = 10)
ggplot(df, aes(Identity, Proportion, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = brewer.pal(11, "RdYlBu")[c(10,2)]) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 32)) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 32)) +
  theme(text = element_text(face = "bold")) +
  theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
  theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
dev.off()

# UMAP velocity stream => on scvelo/scvelo.py

# UMAP velocity length
umap<-read.table("scRNA/diff_ddit3/scvelo/ery_reti_umap_scv.txt")
control.conf<-read.csv("scRNA/diff_ddit3/scvelo/control_confidence_length.csv")
ddit3.conf<-read.csv("scRNA/diff_ddit3/scvelo/ddit3_confidence_length.csv")

control.conf$Condition<-"control"
ddit3.conf$Condition<-"DDIT3"

control.conf$X<-paste0(control.conf$X, "_control")
ddit3.conf$X<-paste0(ddit3.conf$X, "_DDIT3")

conf<-rbind(control.conf, ddit3.conf)

i<-match(colnames(seurat), conf$X)
seurat$velocity_length<-conf[i,3]

pdf("scRNA/diff_ddit3/figures_paper/velocity_length.pdf", useDingbats = F, width = 12, height = 6)
FeaturePlot(seurat, features = "velocity_length", pt.size = 1, cols = rev(brewer.pal(9, "RdYlBu")), split.by = "Condition", order = T) + theme(legend.position="bottom")
dev.off()

# Boxplot velocity length
i<-match(conf$X, colnames(seurat))
conf$CellType<-seurat$CellType[i]

conf<-conf[conf$CellType!="CD34",]
conf<-conf[!is.na(conf$CellType),]

n<-length(unique(conf$CellType))
cc<-data.frame(Stat=rep(NA, n), Pval=rep(NA, n))
ll<-data.frame(Stat=rep(NA, n), Pval=rep(NA, n))

for(n in 1:length(unique(conf$CellType))){
  ix<-unique(conf$CellType)[n]
  df<-conf[conf$CellType==ix,]
  x<-wilcox.test(velocity_confidence~Condition, data = df)
  y<-wilcox.test(velocity_length~Condition, data = df)
  cc[n,1]<-x$statistic
  cc[n,2]<-x$p.value
  ll[n,1]<-y$statistic
  ll[n,2]<-y$p.value
  rownames(cc)[n]<-rownames(ll)[n]<-as.character(ix)
}

cc$Pval.adj<-p.adjust(cc$Pval)
ll$Pval.adj<-p.adjust(ll$Pval)

write.table(cc, file = "scRNA/diff_ddit3/scvelo/confidence_stats.txt", sep = "\t", quote = F)
write.table(ll, file = "scRNA/diff_ddit3/scvelo/length_stats.txt", sep = "\t", quote = F)

pdf("scRNA/diff_ddit3/figures_paper/boxplot_velocity_length_celltype.pdf", width = 12, height = 10)
pp<-ggplot(conf, aes(CellType, velocity_length, fill = Condition)) + geom_boxplot() + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + 
  ggtitle(label = "Length") + labs(y = "Velocity Length")
pp<-pp + theme_bw() + theme(axis.text.x = element_text(angle = 36, hjust = 1)) +
  theme(axis.text = element_text(size = 36, color = "black"), axis.title = element_text(size = 38), plot.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) +
  theme(axis.line.x =  element_line(size = 1.5), axis.ticks.x = element_line(size = 1.5)) + 
  theme(axis.line.y =  element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5))
print(pp)
dev.off()


########################
#####   Figure 5   #####
########################

# GSEA
res<-read.csv(file = "scRNA/diff_ddit3/GSEA_summary_plot_paper.csv", header = T, 
              sep = ";", dec = ",")
groups<-c("Control", "DDIT3")

colnames(res)<-c("GeneSet", "Contrast", "NES", "pval")

sig<-0.4
res$Enriched.group[res$NES>0]<-2
res$Enriched.group[res$NES<0]<-1

res$pval[res$pval>sig]<-sig + 0.01

res$Enriched.group.sig[res$pval>sig]<-0
res$Enriched.group.sig[res$NES>0&res$pval<sig]<-2
res$Enriched.group.sig[res$NES<0&res$pval<sig]<-1

res$Enriched.group<-factor(res$Enriched.group)
res$Enriched.group.sig<-factor(res$Enriched.group.sig)
res$GeneSet<-factor(res$GeneSet, levels = rev(unique(res$GeneSet)))

res$Contrast<-factor(res$Contrast, 
                     levels = unique(res$Contrast))

cols<-brewer.pal(11, "RdYlBu")[c(10,2)]

p<-ggplot(res, aes(Contrast, GeneSet)) + 
  geom_point(aes(fill = Enriched.group.sig, alpha = pval, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES))) +
  scale_shape_identity() 

p<-p + theme_bw() + scale_color_manual(values = cols, labels = groups) +
  scale_fill_manual(values = c("white", cols)) + scale_alpha_continuous(range = c(1,0.001)) +
  scale_size(range = c(2,8)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 14)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 16)) +
  labs(alpha = "P value", size = "NES abs. value", color = "Enriched group") +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
p<-p + theme(text = element_text(face = "plain", colour = "black"))
res.path<-"/home/mainciburu/scRNA/diff_ddit3/figures_paper/gsea_reduced_p0.4.pdf"
pdf(res.path, width = 11.5, height = 5.5, useDingbats = F)
p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
dev.off()

# Pseudotime vs gene trends - GSEA genes and SimiC TFs
gg<-c("SEC61A1", "CBFB", "FAM210B", "ATP5IF1")
tf<-c("SOX6", "ARID4B", "HES6", "TAL1", "MAFG", "E2F2")
gg<-c(gg, tf)
results.control<-readRDS("/home/mainciburu/scRNA/diff_ddit3/palantir_results/gene_trends/control/results_trends_control.rds")
results.ddit3<-readRDS("/home/mainciburu/scRNA/diff_ddit3/palantir_results/gene_trends/ddit3/results_trends_ddit3.rds")

names(results.control)<-"Ery"
names(results.ddit3)<-"Ery"

for(g in gg){
  p<-plot_trends(r1 = results.control, r2 = results.ddit3, gg = g, 
                 bb = "Ery", mode = 211, c1 = "Control", c2 = "DDIT3",
                 col1 = "#4575B4", col2 = "#D73027")
  p <- p + ggtitle(g)
  pdf(paste0("/home/mainciburu/scRNA/diff_ddit3/figures_paper/gene_trends/",
             g, ".pdf"), height = 7, width = 8)
  print(p)
  dev.off()
}

# SimiC heatmap and density plots => analyze_simic_results.r
load('/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_df_auc_KS_D_clust.Rdata')

# Plot AUC distributions per cluster
clusters_2_keep<-c("CFU", "Proerythroblast", "Early_basophilic", "Late_basophilic", 
                   "Polychromatic", "Orthochromatic")

for (clust in clusters_2_keep){
  print(clust)
  plotter <- df_auc[df_auc$cluster_id == clust,]
  pdf(paste0('scRNA/diff_ddit3/figures_paper/AUCs_',clust,'.pdf'), onefile = TRUE, width = 7, height = 7)
  for (tf in unique(df_auc$driver)){
    ptemp<-ggplot(plotter[plotter$driver ==tf,], aes(x=value, fill=.id)) + geom_density(alpha = 0.2, adjust = 1/2, color = "black", size = 0.8) + 
      theme_classic() + scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027")) + 
      theme(title = element_text(size = 40),
            text = element_text(size = 40),
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 40, hjust = 1),
            axis.line = element_line(size = 0.8),
            axis.ticks = element_line(colour = "black", size = 0.8),
            axis.ticks.length=unit(.35, "cm")) +
      labs(fill = "") + 
    ggtitle(tf)
    print(ptemp)
  }
  dev.off()
}

########################
###  Supp Figure 4  ####
########################

# Umap seurat before subseting
seurat.all<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int.rds")
pdf("scRNA/diff_ddit3/figures_paper/umap_beforeSubset_clusters_res.0.4.pdf", 
    useDingbats = F, width = 15, height = 10)
pp<-DimPlot(seurat.all, group.by = "integrated_snn_res.0.4", 
            pt.size = 0.8)
pp<-pp + theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
  theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
print(pp)
dev.off()

# Reticulocytes markers
DefaultAssay(int.ery.reti)<-"RNA"
gg<-c("nCount_RNA", "percent.mito", "percent.ribo", "ATG3", "ATG5", "BNIP3L")
for(g in gg){
  pp<-FeaturePlot(int.ery.reti, features = g, pt.size = 1, cols = c("grey", "red"))

  pp<-pp + theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38)) +
      theme(legend.text = element_text(size = 20)) +
      theme(text = element_text(face = "plain")) + labs(x = "UMAP 1", y = "UMAP 2") +
      theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
      theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
  pdf(paste0("scRNA/diff_ddit3/figures_paper/UMAP_ery_reti_", g, ".pdf"),
      useDingbats = F, width = 15, height = 10)
  print(pp)
  dev.off()
}

# Heatmap singleR => on initial_analysis.r

# VlnPlot DDIT3 (removing reticulocytes)
deg<-readRDS("/home/mainciburu/scRNA/diff_ddit3/deg_condition_ery_clusters_all.rds")
g<-"DDIT3"
Idents(int.ery.reti)<-"CellType"
seurat <- subset(int.ery.reti, idents=c("CD34", "BFU", "CFU", "Proerythroblast", "Early_basophilic", "Late_basophilic", "Polychromatic", "Orthochromatic"))

Idents(seurat)<-"Condition"
seurat$Condition<-factor(seurat$Condition, levels = c("DDIT3", "control"))
df.sig<-deg[deg$gene==g & deg$p_val_adj<0.1,c(5,6)]
df<-data.frame(Identity=levels(seurat$singleR), 
               pval=1,
               significance=factor("NS", levels = c("NS", "*","**")))
i<-na.omit(match(df$Identity, df.sig$cluster))
df.sig<-df.sig[i,]
df$pval[df$Identity%in%df.sig$cluster]<-df.sig$p_val_adj
df$significance[df$pval<0.1]<-"*"
df$significance[df$pval<0.05]<-"**"
pdf(paste0("scRNA/diff_ddit3/figures_paper/vlnplot_", g, ".pdf"), 
    useDingbats = F, width = 10, height = 10)
pp<-VlnPlot(seurat, features = g, group.by = "singleR", split.by = "Condition", assay = "RNA", pt.size = 0) + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + ggtitle(label = g)
n<-nrow(df)
pp<-pp + coord_cartesian(ylim = c(0, 6.5))
ypos<-5.4
pp<- pp + geom_signif(y_position = rep(ypos, n),
                      annotations = df$significance,
                      xmin = 1:n-0.3, xmax = 1:n + 0.3, tip_length = 0.01, 
                      size = 1, textsize = 8)
pp<-pp + theme(axis.text.x = element_text(angle = 36, hjust = 1)) +
  theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38), plot.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) +
  theme(axis.line.x =  element_line(size = 1.5), axis.ticks.x = element_line(size = 1.5)) + 
  theme(axis.line.y =  element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5))
print(pp)
dev.off()



########################
###  Supp Figure 5  ####
########################

# Pseudotime vs gene trends
gg<-c("SAE1", "ZNF574", "WDR18", "VPS11", "MYB", "FOXO3", "HBB", "HBA1", "HBA2", "HBM")

for(g in gg){
  p<-plot_trends(r1 = results.control, r2 = results.ddit3, gg = g, 
                 bb = "Ery", mode = 211, c1 = "Control", c2 = "DDIT3",
                 col1 = "#4575B4", col2 = "#D73027")
  p <- p + ggtitle(g)
  pdf(paste0("/home/mainciburu/scRNA/diff_ddit3/figures_paper/gene_trends/",
             g, ".pdf"), height = 7, width = 8)
  print(p)
  dev.off()
}

########################
### Reviewer answers ###
########################

### UMAP pseudotime
seurat<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds")
pst.control<-read.csv("/home/mainciburu/scRNA/diff_ddit3/palantir_results/control/pseudotime.csv", header = F)
pst.ddit3<-read.csv("/home/mainciburu/scRNA/diff_ddit3/palantir_results/ddit3/pseudotime.csv", header = F)

seurat$plot<-NA
i<-match(colnames(seurat), pst.control$V1)
seurat$plot[!is.na(i)]<-pst.control$V2[i[!is.na(i)]]
i<-match(colnames(seurat), pst.ddit3$V1)
seurat$plot[!is.na(i)]<-pst.ddit3$V2[i[!is.na(i)]]
Idents(seurat)<-"Condition"
cols<-viridis::viridis(30)
pdf(file = "/home/mainciburu/scRNA/diff_ddit3/figures_paper/control_pseudotime.pdf", useDingbats = F,
    width = 6, height = 5)
x<-subset(seurat, idents = "control")
FeaturePlot(x, reduction = "umap", features = "plot", order = T,
            pt.size = 0.2) + ggtitle("Pseudotime") +
  scale_color_gradientn(colours = cols) + labs(colour = "Pseudotime")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

pdf(file = "/home/mainciburu/scRNA/diff_ddit3/figures_paper/ddit3_pseudotime.pdf", useDingbats = F,
    width = 6, height = 5)
x<-subset(seurat, idents = "DDIT3")
FeaturePlot(x, reduction = "umap", features = "plot", order = T,
            pt.size = 0.2) + ggtitle("Pseudotime") +
  scale_color_gradientn(colours = cols) + labs(colour = "Pseudotime")  + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22)) +
  theme(text = element_text(face = "bold")) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

### Gene expression for trend genes
DefaultAssay(seurat)<-"RNA"
gg<-c("SAE1", "ZNF574", "WDR18", "VPS11", "MYB", "FOXO3", "HBB", "HBA1", 
      "HBA2", "HBM", "SEC61A1", "CBFB", "FAM210B", "ATP5IF1",
      "SOX6", "ARID4B", "HES6", "TAL1", "MAFG", "E2F2")

for(g in gg){
  pp<-FeaturePlot(seurat, features = g, pt.size = 1, cols = c("grey", "red"), split.by = "Condition", order = T)

  pp<-pp + theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38)) +
      theme(legend.text = element_text(size = 20)) +
      theme(text = element_text(face = "plain")) + labs(x = "UMAP 1", y = "UMAP 2") +
      theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
      theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
  pdf(paste0("scRNA/diff_ddit3/figures_paper/UMAP_ery_reti_", g, ".pdf"),
      useDingbats = F, width = 12, height = 5)
  print(pp)
  dev.off()
}

### proportion per identity including reticulocytes and myeloid clusters
seurat.all<-readRDS("/home/mainciburu/scRNA/diff_ddit3/seurat_int.rds")
seurat.all$CellType<-as.character(seurat.all$seurat_clusters)
i<-match(colnames(seurat.all), colnames(seurat))
seurat.all$CellType[!is.na(i)]<-as.character(seurat$CellType)[i[!is.na(i)]]
seurat.all$CellType[seurat.all$CellType=="7"]<-"Myeloid_progenitors"
seurat.all$CellType[seurat.all$CellType=="9"]<-"Basophil_progenitors"
seurat.all$CellType<-factor(seurat.all$CellType, levels = c(levels(seurat$CellType), "Myeloid_progenitors", "Basophil_progenitors"))

tt<-prop.table(table(seurat.all$CellType, seurat.all$Condition), margin = 2)
df<-melt(tt)
colnames(df)<-c("Identity", "Condition", "Proportion")
df$Identity<-factor(df$Identity, levels = c(levels(seurat$CellType), "Myeloid_progenitors", "Basophil_progenitors"))
df$Condition<-factor(df$Condition, levels = c("control", "DDIT3"))
pdf("scRNA/diff_ddit3/figures_paper/proportions_myeloid.pdf", 
    useDingbats = F, width = 12, height = 10)
ggplot(df, aes(Identity, Proportion, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = brewer.pal(11, "RdYlBu")[c(10,2)]) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 32)) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 32)) +
  theme(text = element_text(face = "bold")) +
  theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
  theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
dev.off()


### Plot myeloid cluster markers

# Markers
lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")  
gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
granul<-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
mono<-c("LYZ", "CSTA")
baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")
  
markers<-c(lmpp, gmp, granul, mono, dc, baso)
markers<-markers[markers%in%rownames(seurat.all)]
  
seurat.all$CellType<-factor(seurat.all$CellType, levels = c(levels(seurat$CellType), "Myeloid_progenitors", "Basophil_progenitors"))

pdf("scRNA/diff_ddit3/figures_paper/markers_myeloid.pdf", 
    useDingbats = F, width = 10, height = 9)
DotPlot(seurat.all,assay = "RNA",features = rev(markers),dot.min = 0.4,group.by = "CellType",
            col.min = -2,col.max = 2,col=c("#f7e69c","red")) +
            scale_colour_gradient(low = "#f7e69c",high = "red") + theme_classic() + 
            coord_flip() + theme(axis.text.x=element_text(angle=30, hjust = 1),
                                 text=element_text(size=18))
dev.off()

### VlnPlot DDIT3 with myeloid clusters

g<-"DDIT3"
Idents(seurat.all)<-"Condition"
seurat.all$Condition<-factor(seurat.all$Condition, levels = c("DDIT3", "control"))

pdf(paste0("scRNA/diff_ddit3/figures_paper/vlnplot_", g, "_myeloid_clusters.pdf"), 
    useDingbats = F, width = 10, height = 8)
pp<-VlnPlot(seurat.all, features = g, group.by = "CellType", split.by = "Condition", assay = "RNA", pt.size = 0) + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + ggtitle(label = g)
print(pp)
dev.off()

### VlnPlot DDIT3 only myeloid clusters
Idents(seurat.all)<-"CellType"
x<-subset(seurat.all, idents = c("Myeloid_progenitors", "Basophil_progenitors"))

g<-"DDIT3"
Idents(x)<-"Condition"
x$Condition<-factor(x$Condition, levels = c("DDIT3", "control"))

pdf(paste0("scRNA/diff_ddit3/figures_paper/vlnplot_", g, "_ONLY_myeloid_clusters.pdf"), 
    useDingbats = F, width = 10, height = 8)
pp<-VlnPlot(x, features = g, group.by = "CellType", split.by = "Condition", assay = "RNA", pt.size = 0.3) + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + ggtitle(label = g)
print(pp)
dev.off()


### UMAP DDIT3 expression
pp<-FeaturePlot(seurat.all, features = "DDIT3", pt.size = 1, cols = c("grey", "red"), order = T)

pp<-pp + theme(axis.text = element_text(size = 36), axis.title = element_text(size = 38)) +
      theme(legend.text = element_text(size = 20)) +
      theme(text = element_text(face = "plain")) + labs(x = "UMAP 1", y = "UMAP 2") +
      theme(axis.line.x =  element_line(size = 2), axis.ticks.x = element_line(size = 2)) + 
      theme(axis.line.y =  element_line(size = 2), axis.ticks.y = element_line(size = 2)) 
pdf(paste0("scRNA/diff_ddit3/figures_paper/UMAP_", g, "_myeloid_clusters.pdf"),
     useDingbats = F, width = 8, height = 8)
print(pp)
dev.off()


### Velocity confidence and lenght
umap<-read.table("scRNA/diff_ddit3/scvelo/ery_reti_umap_scv.txt")
control.conf<-read.csv("scRNA/diff_ddit3/scvelo/control_confidence_length.csv")
ddit3.conf<-read.csv("scRNA/diff_ddit3/scvelo/ddit3_confidence_length.csv")

control.conf$Condition<-"control"
ddit3.conf$Condition<-"DDIT3"

control.conf$X<-paste0(control.conf$X, "_control")
ddit3.conf$X<-paste0(ddit3.conf$X, "_DDIT3")

conf<-rbind(control.conf, ddit3.conf)

i<-match(colnames(seurat), conf$X)
seurat$velocity_confidence<-conf[i,2]
seurat$velocity_length<-conf[i,3]

pdf("scRNA/diff_ddit3/figures_paper/velocity_confidence.pdf", useDingbats = F, width = 12, height = 6)
FeaturePlot(seurat, features = "velocity_confidence", pt.size = 1, cols = rev(brewer.pal(9, "RdYlBu")), split.by = "Condition", order = T) + theme(legend.position="bottom")
dev.off()
pdf("scRNA/diff_ddit3/figures_paper/velocity_length.pdf", useDingbats = F, width = 12, height = 6)
FeaturePlot(seurat, features = "velocity_length", pt.size = 1, cols = rev(brewer.pal(9, "RdYlBu")), split.by = "Condition", order = T) + theme(legend.position="bottom")
dev.off()


### Boxplots confidence and lenght per cluster
i<-match(conf$X, colnames(seurat))
conf$CellType<-seurat$CellType[i]

conf<-conf[conf$CellType!="CD34",]
conf<-conf[!is.na(conf$CellType),]

n<-length(unique(conf$CellType))
cc<-data.frame(Stat=rep(NA, n), Pval=rep(NA, n))
ll<-data.frame(Stat=rep(NA, n), Pval=rep(NA, n))

for(n in 1:length(unique(conf$CellType))){
  ix<-unique(conf$CellType)[n]
  df<-conf[conf$CellType==ix,]
  x<-wilcox.test(velocity_confidence~Condition, data = df)
  y<-wilcox.test(velocity_length~Condition, data = df)
  cc[n,1]<-x$statistic
  cc[n,2]<-x$p.value
  ll[n,1]<-y$statistic
  ll[n,2]<-y$p.value
  rownames(cc)[n]<-rownames(ll)[n]<-as.character(ix)
}

cc$Pval.adj<-p.adjust(cc$Pval)
ll$Pval.adj<-p.adjust(ll$Pval)

write.table(cc, file = "scRNA/diff_ddit3/scvelo/confidence_stats.txt", sep = "\t", quote = F)
write.table(ll, file = "scRNA/diff_ddit3/scvelo/length_stats.txt", sep = "\t", quote = F)

pdf("scRNA/diff_ddit3/figures_paper/boxplot_velocity_confidence_celltype.pdf", width = 12, height = 10)
pp<-ggplot(conf, aes(CellType, velocity_confidence, fill = Condition)) + geom_boxplot() + ylim(c(0, 1.5)) +
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + 
  ggtitle(label = "Confidence") + labs(y = "Velocity Confidence")
pp<-pp + theme_bw() + theme(axis.text.x = element_text(angle = 36, hjust = 1)) +
  theme(axis.text = element_text(size = 36, color = "black"), axis.title = element_text(size = 38), plot.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) +
  theme(axis.line.x =  element_line(size = 1.5), axis.ticks.x = element_line(size = 1.5)) + 
  theme(axis.line.y =  element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5))
print(pp)
dev.off()


pdf("scRNA/diff_ddit3/figures_paper/boxplot_velocity_length_celltype.pdf", width = 12, height = 10)
pp<-ggplot(conf, aes(CellType, velocity_length, fill = Condition)) + geom_boxplot() + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + 
  ggtitle(label = "Length") + labs(y = "Velocity Length")
pp<-pp + theme_bw() + theme(axis.text.x = element_text(angle = 36, hjust = 1)) +
  theme(axis.text = element_text(size = 36, color = "black"), axis.title = element_text(size = 38), plot.title = element_text(size = 38)) +
  theme(legend.text = element_text(size = 36), legend.title = element_text(size = 38)) +
  theme(text = element_text(face = "plain")) +
  theme(axis.line.x =  element_line(size = 1.5), axis.ticks.x = element_line(size = 1.5)) + 
  theme(axis.line.y =  element_line(size = 1.5), axis.ticks.y = element_line(size = 1.5))
print(pp)
dev.off()

### G2/M S scores

pdf(paste0("scRNA/diff_ddit3/figures_paper/vlnplot_G2MScore_myeloid_clusters.pdf"), 
    useDingbats = F, width = 10, height = 8)
pp<-VlnPlot(seurat.all, features = "G2M.Score", group.by = "CellType", split.by = "Condition", assay = "RNA", pt.size = 0) + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + ggtitle(label = "G2MScore")
print(pp)
dev.off()

pdf(paste0("scRNA/diff_ddit3/figures_paper/vlnplot_SScore_myeloid_clusters.pdf"), 
    useDingbats = F, width = 10, height = 8)
pp<-VlnPlot(seurat.all, features = "S.Score", group.by = "CellType", split.by = "Condition", assay = "RNA", pt.size = 0) + 
  scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027"), guide = guide_legend(reverse = T)) + ggtitle(label = "SScore")
print(pp)
dev.off()





