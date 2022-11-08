setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
task_id = N[task_ID]
dir = "simu_221017_1"
library(tidyverse)
library(igraph)
library(cluster)
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

mu = c(-1,1,0,0)
truth = c(rep(1:2, each = 50), rep(3, 100))
sim_tib = expand_grid(n_feat = c(1000, 10000, 100000),
                      s2n = seq(0.005,0.05,0.005))

source("code/functions/Partition_CIMLR_2.0.R")
kernel_tib = sim_tib %>%
  mutate(tib = map2(n_feat, s2n, function(nf, s2n){
    data = get_data_4clust(n_feat = nf,
                           mu_vec = mu,
                           signal_sd = 1,
                           noise_sd = 1,
                           s2n_ratio = s2n
                           )
    distance = dist2(data, data)
    n = nrow(data)
    k = floor(sqrt(n))
    kernel = kernel_calculation(distance = distance, k = k, sigma = 2)
    diff_kernel = diffusion_enhancement(kernel = kernel, alpha = 0.8, k = k, diffusion_form = "L1")
    cluster_kernel = spectralClustering(kernel, 3)
    cluster_diff_kernel = spectralClustering(diff_kernel, 3)

    cluster_tib = tibble(kernel_type = c("kernel", "diff_kernel"),
                        cluster = list(cluster_kernel, cluster_diff_kernel))

    final_tib = cluster_tib %>%
      mutate(nmi = map_dbl(cluster, function(cluster) compare(cluster, truth, "nmi")),
             adjusted.rand = map_dbl(cluster, function(cluster) compare(cluster, truth, "adjusted.rand")),
             silhouette = map_dbl(cluster, function(cluster){
               sil = silhouette(cluster, distance)
               obj = summary(sil)
               obj$avg.width
             })) %>%
      pivot_longer(cols = c("nmi", "adjusted.rand", "silhouette"),
                   names_to = "metric",
                   values_to = "value"
                   )

    return(final_tib)
  })) %>%
  unnest(tib)

saveRDS(kernel_tib, file = paste(dir, "/save_tib",task_id,".rds", sep = ""))

# kernel_tib %>%
#   ggplot(aes(x = noise_sd, y = nmi, color = kernel_type, group = kernel_type))+
#   geom_line()+
#   geom_point()+
#   facet_grid(n_feat~.)



