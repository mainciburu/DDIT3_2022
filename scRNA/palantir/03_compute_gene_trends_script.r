#########################
# Calculate gene trends 
#########################

#library(feather)
library(hdf5r)
source("/home/mainciburu/scRNA/functions/compute_gene_trends.r")

##### Input data #########
sample<-"ddit3"
sample<-"control"
res_path<-paste0("/home/mainciburu/scRNA/diff_ddit3/palantir_results/", sample, "/")
branch_prob<-read.csv(paste0(res_path, "branch_probs.csv"), row.names = 1)
pseudotime<-read.csv(paste0(res_path, "pseudotime.csv"), row.names = 1)
diff_pot<-read.csv(paste0(res_path, "diff_potential.csv"), header = F, row.names = 1)

f<-H5File$new(paste0(res_path, "imp_df.hdf5"), mode = "r")
g<-f[["imp_df"]][["axis0"]][]
imp_df<-f[["imp_df"]][["block0_values"]][,]
imp_df<-t(imp_df)
colnames(imp_df)<-g
rownames(imp_df)<-rownames(pseudotime)

gene_exprs<-imp_df
lineages<-c("Erythroid_late")
results<-compute_gene_trends(branch_prob = branch_prob, pseudotime = pseudotime,
                             gene_exprs = gene_exprs, lineages = lineages, ncores = 6,
                             res_path = paste0("/home/mainciburu/scRNA/diff_ddit3/palantir_results/gene_trends/", sample, "/"))

saveRDS(results, file = paste0("/home/mainciburu/scRNA/diff_ddit3/palantir_results/gene_trends/", sample, "/results_trends_", sample, "_ery.rds"))
