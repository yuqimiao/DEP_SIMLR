setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
task_id = N[task_ID]
dir = "simu_221012_1"
library(tidyverse)
library(igraph)
# simulation for 4 clusters separates 1/2/34 ----
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

noise_sd_all = c(1, 1.25, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4)
mu = c(-1,1,0,0)
truth = c(rep(1:2, each = 50), rep(3, 100))
sim_tib = expand_grid(n_feat = c(1000, 10000, 100000),
                      noise_sd = seq(1,5,0.25))

source("code/functions/Partition_CIMLR_2.0.R")
kernel_tib = sim_tib %>%
  mutate(kernel_nmi_tib = map2(n_feat, noise_sd, function(nf, sd){
    data = get_data_4clust(n_feat = nf,
                           mu_vec = mu,
                           signal_sd = 1,
                           noise_sd = sd
                           )
    distance = dist2(data, data)
    n = nrow(data)
    k = floor(sqrt(n))
    kernel = kernel_calculation(distance = distance, k = k, sigma = 2)
    diff_kernel = diffusion_enhancement(kernel = kernel, alpha = 0.8, k = k, diffusion_form = "L1")
    cluster_kernel = spectralClustering(kernel, 3)
    cluster_diff_kernel = spectralClustering(diff_kernel, 3)
    nmi_kernel = compare(cluster_kernel, truth, "nmi")
    nmi_diff_kernel = compare(cluster_diff_kernel, truth, "nmi")
    nmi_tib = tibble(kernel_type = c("kernel", "diff_kernel"),
                     nmi = c(nmi_kernel, nmi_diff_kernel))
    return(nmi_tib)
  })) %>%
  unnest(kernel_nmi_tib)

saveRDS(kernel_tib, file = paste(dir, "/save_tib",task_id,".rds", sep = ""))

# kernel_tib %>%
#   ggplot(aes(x = noise_sd, y = nmi, color = kernel_type, group = kernel_type))+
#   geom_line()+
#   geom_point()+
#   facet_grid(n_feat~.)



