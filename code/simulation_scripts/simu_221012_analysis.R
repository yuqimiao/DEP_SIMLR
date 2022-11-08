library(tidyverse)
library(igraph)
library(patchwork)
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

################################################
source("code/functions/visualization_functions.R")

g1 = heatmap_gg(kernel, "Data type 1, kernel")
g2 = heatmap_gg(diff_kernel, "Data type 1, diffused kernel")
g1|g2
ggsave("docs/Extracted_plots/diff_denoise.png")

################################################

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



# explore the data 3 blur effect size ----
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

# explore the eigenspace of the part-cimlr output ----

library(tidyverse)
library(igraph)
# simulation for 4 clusters, data 1 separates 12/3/4, data 2 separates 1/2/34, data 3 separates 1/2/3/4 but with vague division ----
# output
# data matrix with n_sub*n_feat

# simulation par collection

# noise_sd_all = c(1, 1.25, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4)
truth = c(rep(1:2, each = 50), rep(3, 100))

## simulate

n_feat1 = 1000
n_feat2 = 10000
n_feat3 = 100000
noise_sd = 2

mu1 = c(0, 0, 1, -1)
mu2 = c(1,-1,0,0)
mu3 = c(1,2,3,4)


data1 = get_data_4clust(n_feat = n_feat1, mu_vec = mu1, noise_sd = noise_sd)
data2 = get_data_4clust(n_feat = n_feat2, mu_vec = mu2, noise_sd = noise_sd)
data3 = get_data_4clust(n_feat = n_feat2, mu_vec = mu3, noise_sd = noise_sd)
data_list = list(data1, data2, data3)

n_data = length(data_list)

## benchmark
alpha = 1
sigma = 2 # kernel_parameter
diffusion_form = "L1"
n = nrow(data_list[[1]])
k = floor(sqrt(n))
source("code/functions/Partition_CIMLR_2.0.R")
source("code/functions/CIMLR_kernel_input.R")
source("code/functions/visualization_functions.R")
distance_list = lapply(data_list, function(x) dist2(x,x))
kernel_list = lapply(distance_list, function(x) kernel_calculation(distance = x, k = k, sigma = sigma ))
diff_kernel_list = lapply(kernel_list, function(x) diffusion_enhancement(kernel = x, alpha = alpha, k = k , diffusion_form = diffusion_form))

# kernel ----
# CIMLR
res_cimlr = CIMLR_kernel(kernel_list, c = 4, k = k)

# part_CIMLR
res_part_cimlr_obj = part_cimlr(kernel_list, k = k, neig_single = rep(3,n_data), c = 4)
Y = res_part_cimlr_obj$Y
res_part_cimlr = list(S = Y %*% t(Y),
                      cluster = res_part_cimlr_obj$cluster)

kernel_res_list = list(res_cimlr = res_cimlr, res_part_cimlr = res_part_cimlr)

# diff_kernel ----

# CIMLR
res_cimlr_diff = CIMLR_kernel(diff_kernel_list, c = 4, k = k)

# part_CIMLR -- trial 1 neig_single = c(3,3,3)
res_part_cimlr_diff_obj_333 = part_cimlr(diff_kernel_list, k = k, neig_single = c(3,3,3), c = 4)


# part_CIMLR -- trial 2, neig_single = c(3,3,2)
res_part_cimlr_diff_obj_332 = part_cimlr(diff_kernel_list, k = k, neig_single = c(3,3,2), c = 4)

# part_CIMLR -- trial 1 neig_single = c(3,3,4)
res_part_cimlr_diff_obj_334 = part_cimlr(diff_kernel_list, k = k, neig_single = c(3,3,4), c = 4)


# combine to check
Y = res_part_cimlr_diff_obj_334$Y
res_part_cimlr_diff = list(S = Y %*% t(Y),
                           cluster = res_part_cimlr_diff_obj_334$cluster)
diff_kernel_res_list = list(res_cimlr = res_cimlr_diff, res_part_cimlr = res_part_cimlr_diff)

# res_tib
res_tib = expand_grid(kernel = c("kernel", "diff_kernel"),
                      method = c("cimlr", "part_cimlr")) %>%
  mutate(res_list = c(kernel_res_list, diff_kernel_res_list)) %>%
  mutate(nmi = map_dbl(res_list, function(ls) compare(ls$cluster, rep(1:4, each = 50), "nmi")))

res_tib
