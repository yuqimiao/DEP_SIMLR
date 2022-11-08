setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
n = N[task_ID]


source("code/functions/Partition_CIMLR_2.0.R")
source("code/functions/CIMLR_kernel_input.R")
library(SNFtool)
library(tidyverse)

# imputed_ls recipe
data_dir = "data_20220308/pancancer_output"
imputed_ls_dir = "data_20220308/imputed_list"
recipe = tibble(path = list.files(imputed_ls_dir)) %>%
  mutate(cancer = str_sub(path, 1,4)) %>%
  mutate(path = paste(imputed_ls_dir, path, sep = "/"))

path_cur = recipe$path[[n]]
cancer_cur = recipe$cancer[[n]]

# imputed_list read in
imputed_list = readRDS(path_cur)

# parameter setting
n_sample = ncol(imputed_list[[1]])
k = floor(sqrt(n_sample)) # number of neighbors
sigma = 2 # kernel parameter
alpha = 0.8 # diffusion hyper-parameter, explain the contribution of the diffused matrix
diffusion_form = "L1"

kernel_clust_tib_path = paste(data_dir, "/", cancer_cur,"_kernel_clust_tib.rds", sep = "")
if(!file.exists(kernel_clust_tib_path)){
  # distance calculation
  imputed_ls$expr = standardNormalization(t(imputed_ls$expr))
  imputed_ls$me = standardNormalization(t(imputed_ls$me))
  imputed_ls$mu = t(imputed_ls$mu)
  distance_list = lapply(imputed_list, function(x){dist2(x,x)})

  distance_tib = tibble(cancer = cancer_cur,
                        data_type = names(imputed_list),
                        distance = distance_list)
  # kernel generation
  kernel_tib = distance_tib %>%
    mutate(kernel = map(distance_list, function(distance) kernel_calculation(distance, k = k, sigma = sigma)),
           diff_kernel = map(kernel, function(kernel) diffusion_enhancement(kernel, alpha = alpha, k = k, diffusion_form = diffusion_form))) %>%
    pivot_longer(cols = c("kernel","diff_kernel"),
                 names_to = "kernel_type",
                 values_to = "sim_mat")

  # estimate number of factors and single data spectral clustering
  kernel_clust_tib = kernel_tib %>%
    mutate(est_nclust = map_dbl(sim_mat, function(mat) estimateNumberOfClustersGivenGraph(mat, 2:11)[[1]])) %>%
    mutate(cluster = map2(sim_mat, est_nclust, function(S,c) spectralClustering(affinity = S, K = c)))

  saveRDS(kernel_clust_tib, file = kernel_clust_tib_path)
}else{
  kernel_clust_tib = readRDS(kernel_clust_tib_path)
}

# 3 data types integration
final_clust_tib = kernel_clust_tib %>%
  dplyr::select(-distance) %>%
  nest(data = -c(cancer, kernel_type)) %>%
  mutate(tib = map(data, function(dat){
    est_c_single = dat$est_nclust
    sim_list = dat$sim_mat
    integ_tib =  tibble(c = 3:8) %>%
      mutate(SNF = map(c, function(c){
        S_all = SNF(sim_list, k)
        cluster = spectralClustering(S_all,c)
        return(tibble(sim_mat = list(S_all),
                    cluster = list(cluster)))
      })) %>%
      mutate(cimlr = map(c, function(c){
        res = CIMLR_kernel(kernel_list = sim_list,  c = c, k = k)
        return(tibble(sim_mat = list(res$S),
                    cluster = list(res$cluster)))
      })) %>%
      mutate(part_cimlr = map(c, function(c){
        res = part_cimlr(kernel_list = sim_list,k = k, c = c,neig_single = est_c_single, update_neig = F)
        return(tibble(sim_mat  = list(res$Y %*% t(res$Y)),
                    cluster = list(res$cluster)))
      })) %>%
      pivot_longer(cols = c("SNF","cimlr","part_cimlr"),
                   names_to = "method",
                   values_to = "res_tib") %>%
      mutate(data_type = "Integration") %>%
      unnest(res_tib)

    tib = dat %>%
      dplyr::select(data_type, c = est_nclust, sim_mat, cluster) %>%
      mutate(method = "spectral")
    tib = rbind(tib,integ_tib)
    return(tib)
  })) %>%
  unnest(tib)

saveRDS(final_clust_tib, file = paste(data_dir, "/", cancer_cur,"_final_clust_tib.rds", sep = ""))




