setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
task_id = N[task_ID]
dir = "simu_221017_3"
library(tidyverse)
library(igraph)
library(cluster)
# simulation for 4 clusters, data 1 separates 12/3/4, data 2 separates 1/2/34, data 3 separates 1/2/3/4 but with vague division ----
# output
# data matrix with n_sub*n_feat

get_data_4clust = function(n_sub = 200,
                           n_feat = 1000,
                           s2n_ratio = 0.05, # signal to noise ratio
                           mu_vec = c(1, 2, 10, 20), # feature mean for 4 clusters
                           signal_sd = 1,
                           noise_sd = 1){
  data = matrix(0, n_sub, n_feat)
  n_signal = floor(n_feat*s2n_ratio)
  signal_ind = 1:n_signal
  noise_ind = (n_signal+1):n_feat


  n_sub_perc = n_sub/4
  cluster_index = lapply(1:4, function(i) ((i-1)*n_sub_perc+1):(i*n_sub_perc))

  # noise simulation
  data[1:n_sub, noise_ind]=rnorm(n = n_sub*length(noise_ind), mean = 0, sd = noise_sd)


  # signal_simulation
  for(i in 1:4){
    cluster_index_cur = cluster_index[[i]]
    mu_cur = mu_vec[i]
    data[cluster_index_cur, signal_ind]=rnorm(n = n_sub_perc*length(signal_ind), mean = mu_cur, sd = signal_sd)
  }

  return(data)
}


# heatmap_gg(data, "sim1") # heatmap_gg in visualization_functions.R
# simulation par collection

# noise_sd_all = c(1, 1.25, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4)

simu_tib =tibble(scenario = 1:4,
                 n_feat1 = c(1000, 1000, 1000, NA),
                 n_feat2 = c(10000, 10000, NA, 10000),
                 n_feat3 = c(100000, NA, 100000, 100000)) %>%
  mutate(tib = map(scenario, function(s){
    tibble(s2n = seq(0.005,0.05,0.005))
  })) %>% unnest(tib) %>%
  mutate(scenario = seq_along(scenario))


## simulate
scenario_cur = (task_id-1)%/%100+1
simu_tib_cur = simu_tib %>% filter(scenario == scenario_cur)
s2n = simu_tib_cur$s2n
neig_single_vec = c(3,3,4)

n_feat_ls = list(n_feat1 = simu_tib_cur$n_feat1,
                 n_feat2 = simu_tib_cur$n_feat2,
                 n_feat3 = simu_tib_cur$n_feat3)

mu_ls = list(mu1 = c(0, 0, 1, -1),
             mu2 = c(1,-1,0,0),
             mu3 = c(1,2,3,4))

data_list = list()
neig_single = NULL
for(i in 1:3){
  if(!is.na(n_feat_ls[[i]])){
    data = get_data_4clust(n_feat = n_feat_ls[[i]],
                           mu_vec = mu_ls[[i]],
                           s2n_ratio = s2n)
    data_list = c(data_list, list(data))
    neig_single = c(neig_single, neig_single_vec[i])
  }
}


## benchmark
alpha = 1
sigma = 2 # kernel_parameter
diffusion_form = "L1"
n = nrow(data_list[[1]])
k = floor(sqrt(n))
source("code/functions/Partition_CIMLR_2.0.R")
source("code/functions/CIMLR_kernel_input.R")

distance_list = lapply(data_list, function(x) dist2(x,x))
kernel_list = lapply(distance_list, function(x) kernel_calculation(distance = x, k = k, sigma = sigma ))
diff_kernel_list = lapply(kernel_list, function(x) diffusion_enhancement(kernel = x, alpha = alpha, k = k , diffusion_form = diffusion_form))

# kernel ----
# SNF
S_snf = SNF(kernel_list, K = k)
cluster_snf = spectralClustering(S_snf, K = 4)
res_snf = list(S = S_snf,
               cluster = cluster_snf)

# CIMLR
res_cimlr = CIMLR_kernel(kernel_list, c = 4, k = k)

# part_CIMLR
res_part_cimlr = part_cimlr(kernel_list, k = k, neig_single = neig_single, c = 4)
Y = res_part_cimlr$Y
res_part_cimlr = list(S = Y %*% t(Y),
                      cluster = res_part_cimlr$cluster)

# part_CIMLR with updated c
res_part_cimlr_up = part_cimlr(kernel_list, k = k, c = 4, update_neig = T)
Y = res_part_cimlr_up$Y
res_part_cimlr_up = list(S = Y %*% t(Y),
                         cluster = res_part_cimlr_up$cluster)
kernel_res_list = list(res_snf = res_snf, res_cimlr = res_cimlr, res_part_cimlr = res_part_cimlr, res_part_cimlr_up = res_part_cimlr_up)


# diff_kernel ----
# SNF
S_snf = SNF(diff_kernel_list, K = k)
cluster_snf = spectralClustering(S_snf, K = 4)
res_snf = list(S = S_snf,
               cluster = cluster_snf)

# CIMLR
res_cimlr = CIMLR_kernel(diff_kernel_list, c = 4, k = k)

# part_CIMLR
res_part_cimlr = part_cimlr(diff_kernel_list, k = k, neig_single = neig_single, c = 4)
Y = res_part_cimlr$Y
res_part_cimlr = list(S = Y %*% t(Y),
                      cluster = res_part_cimlr$cluster)

# part_CIMLR with updated n_eig
res_part_cimlr_up = part_cimlr(diff_kernel_list, k = k, c = 4, update_neig = T)
Y = res_part_cimlr_up$Y
res_part_cimlr_up = list(S = Y %*% t(Y),
                         cluster = res_part_cimlr_up$cluster)
diff_kernel_res_list = list(res_snf = res_snf, res_cimlr = res_cimlr, res_part_cimlr = res_part_cimlr, res_part_cimlr_up = res_part_cimlr_up)

# res_tib

ave_dist = Reduce("+",distance_list)/length(distance_list)
res_tib = expand_grid(kernel = c("kernel", "diff_kernel"),
                      method = c("snf", "cimlr", "part_cimlr","part_cimlr_up")) %>%
  mutate(res_list = c(kernel_res_list, diff_kernel_res_list)) %>%
  mutate(nmi = map_dbl(res_list, function(ls) compare(ls$cluster, rep(1:4, each = 50), "nmi")),
         adjusted.rand = map_dbl(res_list, function(ls) compare(ls$cluster, rep(1:4, each = 50), "adjusted.rand")),
         silhouette = map_dbl(res_list, function(ls){
           sil = silhouette(ls$cluster, ave_dist)
           obj = summary(sil)
           obj$avg.width
         })
  ) %>%
  pivot_longer(cols = c("nmi", "adjusted.rand", "silhouette"),
               names_to = "metric",
               values_to = "value"
  )

save_tib = simu_tib_cur %>%
  mutate(res_tib = list(res_tib)) %>%
  unnest(res_tib) %>%
  dplyr::select(-res_list)

saveRDS(save_tib, file = paste(dir, "/save_tib",task_id,".rds", sep = ""))






















