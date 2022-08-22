
library(reticulate)
reticulate::use_python("/home/mainciburu/scRNA/simic/myenv/bin/python3")
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggridges)
library(stringr)
library(xlsx)


############## input ####################
data_root_dir <-"/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/"
analysis_root_dir <- "/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/"
plot_root_dir <-"/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/plots/"
cluster_names <- c('control', 'DDIT3')
assignment<-as.character(seq(0,length(cluster_names)-1))
names(assignment) <- cluster_names
plasma <- viridis(50, direction = 1, option = "C")

# load weights and AUC
Ws_f <- paste0(analysis_root_dir, "ddit3_weights_filtered_BIC.pickle")
AUC_files<-dir(path = analysis_root_dir, pattern = "BIS.csv", full.names = F)

weights_ <- py_load_object(filename =Ws_f)

r<-0
AUCs<-list()
for(AUC_f in AUC_files){
  r<-r+1
  a <- read.table(paste0(analysis_root_dir, AUC_f), header=T, sep='\t')
  rownames(a) <- a$X
  a <- a[,!colnames(a) %in% c('X')]
  AUCs[[assignment[r]]]<-a
}

saveRDS(AUCs, '/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_AUCs.rds')

# extract weight matrices (100TF x 1000 targets)
Ws<-list()
for(i in seq(1:length(weights_$weight_dic))){
  W<-weights_$weight_dic[[i]][-nrow(weights_$weight_dic[[i]]),]
  #W<-scale(W)
  rownames(W)<-weights_$TF_ids
  colnames(W)<-weights_$query_targets
  Ws[[cluster_names[i]]]<-W
}

saveRDS(Ws, '/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_Ws.rds')

# Load seurat data
seurat_data <- readRDS('/home/mainciburu/scRNA/diff_ddit3/seurat_int_ery_reti.rds')

# Load AUCs and Weights as .rds
AUCs<-readRDS('/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_AUCs.rds')
Ws<-readRDS('/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_Ws.rds')

#######################

############## Weight plots ################################
# Histogram of R squared per target
# There should be > 700 targets and mean R2 > 0.7 - 0.8
unselected_targets <- list()
unselected_targets[[assignment[['control']]]] <- weights_$query_targets[which(weights_$adjusted_r_squared[[assignment[['control']]]] < 0.7)]
unselected_targets[[assignment[['DDIT3']]]] <- c(weights_$query_targets[which(weights_$adjusted_r_squared[[assignment[['DDIT3']]]] < 0.7)])

pdf(paste0(plot_root_dir, 'Ws_hist_filtered_BIC.pdf'))
selectedTargets <- which(weights_$adjusted_r_squared[[assignment[['control']]]] > 0.7)
hist(weights_$adjusted_r_squared[[assignment[['control']]]], col='grey', breaks=100, main = paste0('Targets selected:: ', length(selectedTargets), ', mean:: ', mean(weights_$adjusted_r_squared[[assignment[['control']]]][selectedTargets])))
selectedTargets <- which(weights_$adjusted_r_squared[[assignment[['DDIT3']]]] > 0.7)
hist(weights_$adjusted_r_squared[[assignment[['DDIT3']]]], col='grey', breaks=100, main = paste0('Targets selected:: ', length(selectedTargets), ', mean:: ', mean(weights_$adjusted_r_squared[[assignment[['DDIT3']]]][selectedTargets])))
dev.off()

# weight matrix to df with weight per TF and target
df <- ldply (Ws, data.frame)
df$driver<-unlist(lapply(Ws,function(x)rownames(x)))
df_w <- melt(df ,variable.name = 'target')

# Weights barplot per TF and condition
pdf_name_ws <- paste0(plot_root_dir, "w_filtered_BIC.pdf")
TF_2_remove <- c()
pdf(pdf_name_ws, onefile = TRUE, width=20)
for(drv in unique(df_w$driver)){
  cat(paste0(drv, '\n'))
  # get the target genes for each TF
  tmp_plotter <- df_w[df_w$driver == drv,]       # target weights for the TF
  # tmp_plotter$value <- scale(tmp_plotter$value, center=FALSE)
  # remove the genes tf with less than 0.7 R2
  tmp_plotter_control <- tmp_plotter[tmp_plotter$.id == 'control' & !tmp_plotter$target %in% c(unselected_targets[[assignment[['control']]]]) ,]
  tmp_plotter_DDIT3 <- tmp_plotter[tmp_plotter$.id == 'DDIT3' & !tmp_plotter$target %in% c(unselected_targets[[assignment[['DDIT3']]]]) ,]
  tmp_plotter <- rbind(tmp_plotter_control, tmp_plotter_DDIT3)
  #
  tmp_plotter <- tmp_plotter[order(abs(tmp_plotter$value), decreasing=T),]    # sort by decreasing weight
  bests <- unique(tmp_plotter$target)[1:100]    # choose 100 top targets
  tmp_plotter <- tmp_plotter[tmp_plotter$target %in% bests,]
  tmp_plotter$target <- factor(tmp_plotter$target, levels = unique(tmp_plotter$target))
  name <- ifelse(max(abs(tmp_plotter$value)) > 0.001, drv, paste0(drv, '*'))    # remove TF if maximum target weight is < 0.001
  if(grepl('\\*', name)){TF_2_remove <- c(TF_2_remove, drv)}
  p <- ggplot(tmp_plotter, aes(x=target, y=value, fill=.id,  palette = "jco")) + 
    geom_bar(stat='identity', position='dodge', color='black') + 
    scale_fill_manual(values = c('darkorange', 'dodgerblue')) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(name) 
  print(p)
}
dev.off()

# Weight barplot per target and condition
pdf(paste0(plot_root_dir, 'Target_w_filtered_BIC.pdf'), onefile = TRUE, width=20)
plot_counter <- 1
plot_list <- list()
for(tgt in sort(as.character(unique(df_w$target)))){
  cat(paste0(tgt, '\n'))
  # get the target genes for each TF
  if(tgt %in% intersect(unselected_targets[[assignment[['control']]]], unselected_targets[[assignment[['DDIT3']]]])){
    next
  }else if(tgt %in% c(unselected_targets[[assignment[['control']]]], unselected_targets[[assignment[['DDIT3']]]])){
    
    ifelse( tgt %in% unselected_targets[[assignment[['control']]]] , tmp_plotter <- df_w[df_w$target == tgt & df_w$.id == 'DDIT3',] , tmp_plotter <- df_w[df_w$target == tgt & df_w$.id == 'control',])
  }else{
    tmp_plotter <- df_w[df_w$target == tgt,]
  }
  # tmp_plotter$value <- scale(tmp_plotter$value, center=FALSE)
  tmp_plotter <- tmp_plotter[order(abs(tmp_plotter$value), decreasing=T),]
  assign(paste0('p', plot_counter), ggplot(tmp_plotter, aes(x=reorder(driver, -abs(value)), y=value, fill=.id,  palette = "jco")) + 
           geom_bar(stat='identity', position='dodge', color='black') + 
           scale_fill_manual(values = c('darkorange', 'dodgerblue')) + 
           theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(tgt) )
  if(plot_counter == 4){
    grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
    plot_list <- list()
    plot_counter <- 1
    print('plotted')
  }else{
    plot_counter <- plot_counter +1
  }
}
dev.off()

############################################################################################################

############## AUC plots ################################

# Separate AUC by state
Idents(seurat_data)<-seurat_data$Condition
AUCs_by_state<-list()
state_specific_AUCs<-NULL
for(i in 1:length(cluster_names)){
  cell_names_i<-colnames(subset(seurat_data, idents = cluster_names[i]))
  AUCs_cell_i<-AUCs[[assignment[i]]]
  AUCs_by_state[[cluster_names[i]]]<-AUCs_cell_i[cell_names_i,]
  state_specific_AUCs<-rbind(state_specific_AUCs, AUCs_cell_i[cell_names_i,])
}

# !!!I had remove 6 and 7 cells from control and DDIT3. They are NA in these AUC matrices (not present in original AUCs)

# Add cluster information
df <- ldply(AUCs_by_state, data.frame)
df$cell_id<-rownames(state_specific_AUCs)
clusters.df<-tibble::rownames_to_column(as.data.frame(seurat_data$singleR),"cell_id")
clusters.df<-dplyr::rename(clusters.df,cluster_id='seurat_data$singleR')
df_w_cluster<-merge(df,clusters.df,by="cell_id")
df_auc <- melt(df_w_cluster ,  id.vars = c('.id','cell_id', 'cluster_id'), variable.name = 'driver', stringsAsFactors =F)
# ! NAs removed

# Plot cell number per cluster
tmp <- with(df_w_cluster, table(cluster_id, .id))
tmp <- melt(tmp)
pdf(paste0(plot_root_dir,'cluster_pop.pdf'))
ggplot(tmp, aes(x=cluster_id, y=value, fill=.id)) + 
  geom_bar(stat='identity' , position='dodge',  color='black') + 
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_manual(values = c('darkorange', 'dodgerblue')) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
dev.off()

# clusters with cells in both conditions
clusters_2_keep <- as.data.frame.matrix(with(df_auc, table(cluster_id, .id)))
clusters_2_keep$cluster_id <- rownames(clusters_2_keep)
clusters_2_keep <- clusters_2_keep[clusters_2_keep$control >0 & clusters_2_keep$DDIT3 >0, 'cluster_id']

# Separate AUCs by condition and cluster
AUCs_by_CellType_control<-list()
state_specific_AUCs_CellType_control<-NULL
AUCs_by_CellType_DDIT3<-list()
state_specific_AUCs_CellType_DDIT3<-NULL
for(clust in unique(clusters_2_keep)){
  Idents(seurat_data) <- seurat_data$singleR
  cell_names_i<-colnames(subset(seurat_data, idents = clust))
  # control
  AUCs_cell_i<-AUCs[['0']][cell_names_i,]
  AUCs_by_CellType_control[[clust]]<-AUCs_cell_i
  state_specific_AUCs_CellType_control<-rbind(state_specific_AUCs_CellType_control, AUCs_cell_i[cell_names_i,])
  # DDIT3
  AUCs_cell_i<-AUCs[['1']][cell_names_i,]
  AUCs_by_CellType_DDIT3[[clust]]<-AUCs_cell_i
  state_specific_AUCs_CellType_DDIT3<-rbind(state_specific_AUCs_CellType_DDIT3, AUCs_cell_i[cell_names_i,])
}

df_auc_common <- df_auc[df_auc$cluster_id %in% clusters_2_keep, ]

# Test AUC distribution per TF among cells in the same cluster between conditions
KS_D_clust<-NULL
for(i in unique(df_auc_common$cluster_id)){    # for each cluster
  #browser()
  wauc.ci.df <- dplyr::filter(df_auc_common,cluster_id==i)
  KS_D<-NULL
  for (tf in unique(df_auc_common$driver)){    # for each TF
    wauc.ci.tf.soldier <- dplyr::filter(wauc.ci.df,.id=='control'&driver==tf)$value
    wauc.ci.tf.forager <- dplyr::filter(wauc.ci.df,.id=='DDIT3'&driver==tf)$value
    KS_D<-append(KS_D,ks.test(wauc.ci.tf.soldier,ecdf(wauc.ci.tf.forager))$statistic)
  }
  KS_D_clust<-cbind(KS_D_clust,KS_D)
}
colnames(KS_D_clust)<-unique(df_auc_common$cluster_id)
rownames(KS_D_clust)<-unique(df_auc_common$driver)

# save df_auc and KS_D_clust for plotting
save(df_auc, KS_D_clust, file = "/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/ddit3_filtered_BIC_df_auc_KS_D_clust.Rdata")

# Plot AUC distributions per cluster
for (clust in clusters_2_keep){
  print(clust)
  plotter <- df_auc[df_auc$cluster_id == clust,]
  pdf(paste0(plot_root_dir, 'AUCs_Cluster__',clust,'__Filtered_BIC.pdf'), onefile = TRUE)
  plot_counter <- 1
  for (tf in unique(plotter$driver)){
    ptemp<-ggplot(plotter[plotter$driver ==tf,], aes(x=value, fill=.id)) + geom_density(alpha = 0.2, adjust = 1/2) + 
            theme_classic() + scale_fill_manual(values = c(control="#4575B4", DDIT3 = "#D73027")) + 
            theme(title = element_text(size = 16),
                  text = element_text(size = 14))
            ggtitle(paste0(tf, '    ', KS_D_clust[tf, clust])) 
    assign( paste0('p', plot_counter), ptemp)
    if(plot_counter == 4){
      grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
      plot_counter <- 1
      print('plotted')
    }else{
      plot_counter <- plot_counter +1
    }
  }
  grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
  dev.off()
}

# Plot heatmap with KS_D values per cluster
# TFs to remove???
separable_TFs <- which(apply(KS_D_clust,1,mean)>0.8)
KS_D_clust_filtered<-KS_D_clust[-separable_TFs,]

clust_order_asc<-names(sort(apply(KS_D_clust, 2, mean)))    # order clusters by KS_D
clust_order<-levels(seurat_data$singleR)[-(1:2)]
KS_D_clust<-KS_D_clust[,clust_order]

KS_D_df<-tibble::rownames_to_column(as.data.frame(KS_D_clust), var = 'driver')
KS_D_df<-melt(KS_D_df,variable.name = "cluster_id")

library(pheatmap)
library(gridExtra)
library(stringr)
# !! Go to lines 68 - 110 in TF_2_remove is not created
KS_D_clust <- KS_D_clust[!rownames(KS_D_clust) %in% TF_2_remove,]
pdf(paste0(plot_root_dir, 'Umap_KsAuc_HeatMap_Filtered_BIC.pdf'), onefile = TRUE)
pheatmap(KS_D_clust, color=plasma, cluster_cols =F, fontsize=5, angle_col =45)
pheatmap(KS_D_clust, color=plasma, cluster_rows =F, cluster_cols =F, fontsize=5, angle_col =45)
pheatmap(KS_D_clust, color=plasma, cluster_rows =F, fontsize=5, angle_col =45)
pheatmap(KS_D_clust, color=plasma, fontsize=5, angle_col =45)
dev.off()

# Plot KS_D value distribution per cluster
KS_D_df <- KS_D_df[!KS_D_df$driver %in% TF_2_remove,]
pdf(paste0(plot_root_dir, 'Umap_KsAuc_Ridges_Filtered_BIC.pdf'), onefile = TRUE)
ggplot(KS_D_df, aes(x = value, y = cluster_id, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "KS Dist, F vs S", option = "C") + 
  labs(title = 'Distribution of KS distance between distr of wAUC of control and DDIT3.')
dev.off()

# Plot KS_D value binarized distribution per cluster
pdf(paste0(plot_root_dir, "Umap_KsAuc_BinRidges_Filtered_BIC.pdf"), onefile = TRUE)
ggplot(KS_D_df, aes(x = value, y = cluster_id, fill = stat(x))) +
  geom_density_ridges_gradient(stat = "binline", binwidth= 0.02, draw_baseline = FALSE, scale = 3, rel_min_height = 0.01) + 
  xlim(0,1) + scale_fill_viridis_c(name = "KS Dist, L vs H", option = "C") + 
  labs(title = 'Distribution of KS distance between distr of wAUC of control and DDIT3.')
dev.off()

# Plot KS_D value distribution per cluster (with violin plot)
pdf(paste0(plot_root_dir, "Umap_KsAuc_Boxplots_Filtered_BIC.pdf"), onefile = TRUE)
ggplot(KS_D_df, aes(x = cluster_id, y = value, fill = stat(x))) +
  geom_violin() + scale_fill_viridis_c(name = "KS Dist, L vs H", option = "C") + 
  labs(title = 'Distribution of KS distance between distr of wAUC of control and DDIT3.')
dev.off()

# Plot AUC distribution for every cluster
pdf_name_auc <- paste0(plot_root_dir, "auc_Filtered_BIC.pdf")
pdf(pdf_name_auc, onefile = TRUE)
plot_counter <- 1
for (drv in unique(df_auc$driver)){
  print(drv)
  assign(paste0('p', plot_counter), ggplot(df_auc[df_auc$driver == drv,] , aes(x=value, fill=.id)) + geom_density(alpha = 0.2) + theme_classic() + ggtitle(drv))
  if(plot_counter == 2){
    grid.arrange(p1,p2, ncol=2)
    plot_counter <- 0
  }
  plot_counter <- plot_counter +1
}
dev.off()








