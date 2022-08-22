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





data_folder <- '/home/mainciburu/scRNA/simic/l1_0.01_l2_0.1/results/'


files <- c('ddit3_weights.pickle')
for (fl in files){
    ifelse(file.exists(paste0(data_folder, fl)), print('ok'), print(fl))
}



for(fl in files){
    #name <- str_extract(fl, '^[\\w]+(?=_L10)')
    name<-fl
    pdf(paste0(data_folder, name, '_Weigths_BIC.pdf'))
    weigths <- py_load_object(paste0(data_folder, fl)) # load weigths dictionary
    print(fl)
    for (name_w in names(weigths[['weight_dic']])){     # for each condition (control / DDIT3)
        # work with weight matrix for one condition, with (100TF + 1 misterious row) x 1000 targets
        print(name_w)
        max_l <- c()
        bottom <- weigths$weight_dic[[name_w]][101,]    
        tmp <- as.data.frame(weigths$weight_dic[[name_w]][-nrow(weigths$weight_dic[[name_w]]),])   # remove row 101
        colnames(tmp) <- weigths[['query_targets']]
        rownames(tmp) <- weigths[['TF_ids']]
        tmp <- scale(tmp, center=FALSE)
        for( target in weigths[['query_targets']]){
            all_tf <- rownames(tmp)[order(abs(tmp[,target]), decreasing=T)]  # order TFs from greater to smaller association with that target
            l <- 1
            tf_2_keep <- all_tf[1:l]     # Keep the 1st TF
            while ( sum(tmp[tf_2_keep, target]^2) /sum(tmp[,target]^2)  < 0.9 ){     # Add TFs until this threshold
                l <- l +1
                tf_2_keep <- all_tf[1:l]
            }
            max_l <- c(max_l, l)     # Store TF number for that target
            tmp[!rownames(tmp) %in% tf_2_keep, target] <- 0       # Change to 0 weights for that target and TFs not included in tf_2_keep
            
        }
        tmp <- rbind(tmp, bottom)       # Add again row 101
        tmp <- as.matrix(tmp)
        weigths$weight_dic[[name_w]] <- tmp
        hist(max_l, breaks = 20, main = paste0(name, '_Weights_',name_w ), col = 'darkgrey')     # histogram of TF number per target
    }
    new_fl_name = paste0(data_folder,sub('.pickle', '_filtered_BIC.pickle', fl))
    reticulate::py_save_object(weigths, filename = new_fl_name)
    dev.off()
    print(new_fl_name)
}
