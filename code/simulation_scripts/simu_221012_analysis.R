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

source("code/functions/Partition_CIMLR_2.0.R")

# effect of diffusion ----
# case 1 ----
data = get_data_4clust(mu_vec = c(-1,1,0,0),
                       noise_sd = 1)
truth = c(rep(1,50), rep(2, 50), rep(3, 100))
n = nrow(data)
k = floor(sqrt(n))
distance = dist2(data, data)
kernel = kernel_calculation(distance,k = k, sigma = 2)
diff_kernel = diffusion_enhancement(kernel, alpha = 0.8, k = k,diffusion_form = "L1")

cluster_kernel = spectralClustering(kernel,3)
cluster_diff_kernel = spectralClustering(diff_kernel,3)

compare(cluster_kernel, truth, "nmi")
compare(cluster_diff_kernel, truth, "nmi")

# case 2, reduce noise ----
data = get_data_4clust(mu_vec = c(-1,1,0,0),
                       noise_sd = 0.5)
truth = c(rep(1,50), rep(2, 50), rep(3, 100))
n = nrow(data)
k = floor(sqrt(n))
distance = dist2(data, data)
kernel = kernel_calculation(distance,k = k, sigma = 2)
diff_kernel = diffusion_enhancement(kernel, alpha = 0.8, k = k,diffusion_form = "L1")

cluster_kernel = spectralClustering(kernel,3)
cluster_diff_kernel = spectralClustering(diff_kernel,3)

compare(cluster_kernel, truth, "nmi")
compare(cluster_diff_kernel, truth, "nmi")

# case 3,increase dimension ----
data = get_data_4clust(mu_vec = c(-1,1,0,0),
                       noise_sd = 1,
                       n_feat = 10000)
truth = c(rep(1,50), rep(2, 50), rep(3, 100))
n = nrow(data)
k = floor(sqrt(n))
distance = dist2(data, data)
kernel = kernel_calculation(distance,k = k, sigma = 2)
diff_kernel = diffusion_enhancement(kernel, alpha = 0.8, k = k,diffusion_form = "L1")

cluster_kernel = spectralClustering(kernel,3)
cluster_diff_kernel = spectralClustering(diff_kernel,3)

compare(cluster_kernel, truth, "nmi")
compare(cluster_diff_kernel, truth, "nmi")



# explore the data effect size ----
source("code/functions/visualization_functions.R")
data = get_data_4clust(mu_vec = c(1,2,3,4),
                       noise_sd = 1)
truth = c(rep(1,50), rep(2, 50), rep(3, 50), rep(4, 50))
n = nrow(data)
k = floor(sqrt(n))
distance = dist2(data, data)
kernel = kernel_calculation(distance,k = k, sigma = 2)
diff_kernel = diffusion_enhancement(kernel, alpha = 0.8, k = k,diffusion_form = "L2")

heatmap_gg(eigen(kernel)$vector[,1:10],"kernel_vec")

s0 = kernel
diag(s0) = mean(s0)
heatmap_gg(s0, "4 clust blur")

s0 = diff_kernel
diag(s0) = mean(s0)
heatmap_gg(s0, "4 clust blur, diffusion")

cluster_kernel = spectralClustering(kernel,4)
cluster_diff_kernel = spectralClustering(diff_kernel,4)

compare(cluster_kernel, truth, "nmi")
compare(cluster_diff_kernel, truth, "nmi")

