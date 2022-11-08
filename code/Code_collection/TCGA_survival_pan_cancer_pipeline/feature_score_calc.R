setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
n = N[task_ID]
library(SNFtool)
library(tidyverse)
source("code/functions/feature_score.R")

# dir = "/Volumes/sw2206/yuqi"
dir = "."

cluster_tib_dir = paste(dir, "/","data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer", sep = "")

recipe = tibble(cancer = c("KIRC", "KIRP","LUAD"),
                cluster_tib_path = c(paste(dir, "/", "data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/KIRC_final_clust_tib.rds", sep = ""),
                                     paste(dir, "/", "data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/KIRP_final_clust_tib.rds", sep = ""),
                                     paste(dir, "/", "data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/LUAD_final_clust_tib.rds", sep = "")),

                imputed_ls_path = c(paste(dir, "/","data_20220308/imputed_list/KIRC_imputDat_ls.rds",sep = ""),
                                    paste(dir, "/","data_20220308/imputed_list/KIRP_imputDat_ls.rds",sep = ""),
                                    paste(dir, "/","data_20220308/imputed_list/LUAD_imputDat_ls.rds",sep = "")),
                c = c(3,4,4))


imputed_ls_path = recipe$imputed_ls_path[[n]]
cluster_tib_path = recipe$cluster_tib_path[[n]]
c_cur = recipe$c[[n]]
cancer_cur = recipe$cancer[[n]]

imputed_ls = readRDS(imputed_ls_path)
cluster_tib = readRDS(cluster_tib_path)

W = (cluster_tib %>%
       filter(c == c_cur) %>%
       filter(method == "part_cimlr") %>%
       filter(kernel_type == "diff_kernel") %>%
       pull(sim_mat))[[1]]

feat_rank_ls = lapply(imputed_ls, function(X){
  feat_rank = SIMLR_Feature_Ranking(A = W, X = X)
  return(feat_rank)
})



saveRDS(feat_rank_ls, file = paste(cluster_tib_dir,"/",cancer_cur,"_feature_rank.rds", sep = ""))
