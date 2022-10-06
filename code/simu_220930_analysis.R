# mds plot for simulation result
cur_id = 304
sim_res = tib_all %>%
  filter(noise_sd == 2) %>%
  filter(simu_id == cur_id) %>%
  unnest(res_tib)

part_S_kernel = (sim_res %>% filter(kernel == "kernel" & method == "part_cimlr") %>% pull(res_list))[[1]]$S
part_S_diff_kernel = (sim_res %>% filter(kernel == "diff_kernel" & method == "part_cimlr") %>% pull(res_list))[[1]]$S

part_S_kernel_cl = (sim_res %>% filter(kernel == "kernel" & method == "part_cimlr") %>% pull(res_list))[[1]]$cluster
part_S_diff_kernel_cl = (sim_res %>% filter(kernel == "diff_kernel" & method == "part_cimlr") %>% pull(res_list))[[1]]$cluster

source("code/functions/Partition_CIMLR_2.0.R")
part_S_kernel_dist = dist_kernels(part_S_kernel)
part_S_diff_kernel_dist = dist_kernels(part_S_diff_kernel)

source("code/functions/visualization_functions.R")
# Labeled by the cluster result
mds_kernel = as_tibble(cmdscale(part_S_kernel_dist, 2))%>%
  mutate(cluster = part_S_kernel_cl) %>%
  ggplot(aes(x = V1, y = V2, color = factor(cluster))) +
  geom_point()+
  labs(title = "MDS plot of dist(YY^T) from part-cimlr-kernel ",
       color = "cluster label")

mds_diff_kernel = as_tibble(cmdscale(part_S_diff_kernel_dist, 2))%>%
  mutate(cluster = part_S_diff_kernel_cl) %>%
  ggplot(aes(x = V1, y = V2, color = factor(cluster))) +
  geom_point()+
  labs(title = "MDS plot of dist(YY^T) from part-cimlr-diff_kernel ",
       color = "cluster label")

# Labeled by the cluster truth
mds_kernel_t = as_tibble(cmdscale(part_S_kernel_dist, 2))%>%
  mutate(cluster = rep(1:4, each = 50)) %>%
  ggplot(aes(x = V1, y = V2, color = factor(cluster))) +
  geom_point()+
  labs(title = "MDS plot of dist(YY^T) from part-cimlr-kernel ",
       color = "truth")

mds_diff_kernel_t = as_tibble(cmdscale(part_S_diff_kernel_dist, 2))%>%
  mutate(cluster = rep(1:4, each = 50)) %>%
  ggplot(aes(x = V1, y = V2, color = factor(cluster))) +
  geom_point()+
  labs(title = "MDS plot of dist(YY^T) from part-cimlr-diff_kernel ",
       color = "truth")

(mds_kernel|mds_diff_kernel)/
(mds_kernel_t|mds_diff_kernel_t)

# check the kernels and diff_kernels ----
library(tidyverse)
library(igraph)
# simulation for 4 clusters, data 1 separates 1/2A/2B, data 2 separates 1A/1B/2 ----
# output:
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

# get kernel and diff kernel ----


## simulate
n_feat1 = 1000
n_feat2 = 100000
mu1 = c(0, 0, 1, -1)
mu2 = c(1,-1,0,0)
noise_sd = 3.25

data1 = get_data_4clust(n_feat = n_feat1, mu_vec = mu1, noise_sd = noise_sd)
data2 = get_data_4clust(n_feat = n_feat2, mu_vec = mu2, noise_sd = noise_sd)

data_list = list(data1, data2)

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


S1 = kernel_list[[1]]
S2 = kernel_list[[2]]
diag(S1) = mean(S1)
diag(S2) = mean(S2)
heatmap_gg(S1,"data kernel")
heatmap_gg(S2,"data2, kernel")
heatmap_gg(S1+S2,"kernel sum")

heatmap_gg(eigen(S1)$vectors[,1:5],"S1 first 5")
heatmap_gg(eigen(S2)$vectors[,1:5],"S2 first 5")

est_nclust(S1)
est_nclust(S2)



L1 = normalized_GL(S1)
plot(2:10, sort(eigen(L1)$values)[2:10], main = "GL 1 1-5 smallest eigenval")

L2 = normalized_GL(S2)
plot(2:10, sort(eigen(L2)$values)[2:10], main = "GL 2 1-5 smallest eigenval")

heatmap_gg(eigen(L1)$vectors[,(n-5):n],"GL 1-5 smallest eigenvec")
heatmap_gg(eigen(L2)$vectors[,(n-5):n],"GL 1-5 smallest eigenvec")

compare(kmeans(eigen(L1)$vectors[,(n-2):n],3, nstart = 200)$cluster, c(rep(1,100), rep(2:3, each = 50)), "nmi")
compare(kmeans(eigen(L1)$vectors[,(n-1):n],3, nstart = 200)$cluster, c(rep(1,100), rep(2:3, each = 50)), "nmi")

compare(kmeans(eigen(L2)$vectors[,(n-2):n],3, nstart = 200)$cluster, c(rep(1:2, each = 50), rep(3,100)), "nmi")
compare(kmeans(eigen(L2)$vectors[,(n-1):n],3, nstart = 200)$cluster, c(rep(1:2, each = 50), rep(3,100)), "nmi")

# cbind Y
Y = cbind(eigen(L1)$vectors[,(n-1):n], eigen(L2)$vectors[,(n-1):n])
compare(kmeans(Y, 4, nstart = 200)$cluster, rep(1:4, each = 50), "nmi")


Y = cbind(eigen(L1)$vectors[,(n-2):n], eigen(L2)$vectors[,(n-2):n])
compare(kmeans(Y, 4, nstart = 200)$cluster, rep(1:4, each = 50), "nmi")

# learned Y, no weight
c_single = c(2,2)

Fs = list(data1 = eigen(L1)$vectors[,(n-c_single[1]+1):n], data2 = eigen(L2)$vectors[,(n-c_single[2]+1):n])
Ls = lapply(Fs, function(f) diag(1,n)-f %*% t(f) *2)
Ly = Reduce("+", Ls)

c = 4
Y = eigen(Ly)$vectors[, (n-c+1):n]
plot(2:10, sort(eigen(Ly)$values)[2:10], main = "eigvec from learned Y")
heatmap_gg(Y, "learned Y, no weight")
compare(kmeans(Y, 4, nstart = 200)$cluster, rep(1:4, each = 50), "nmi")

c = 3
Y = eigen(Ly)$vectors[, (n-c+1):n]
compare(kmeans(Y, 4, nstart = 200)$cluster, rep(1:4, each = 50), "nmi")

# update diffusion_kernel makes worse performance?
kernel_list = kernel_list
c=4
update_c = T
num_eig = 2:10
rho = 1

# kernel distance calculation ----
S = length(kernel_list) # number of data types
n = nrow(kernel_list[[1]]) # number of samples
kernel_distance_list = lapply(kernel_list, function(kernel) dist_kernels(kernel))

# initialization ----
# par initial
initial_list = lapply(kernel_distance_list, function(distance) initial(distance, k = k))
Z0 = lapply(kernel_distance_list, function(x){
  x = as.matrix(x)
  max(x)-x
})

# estimate eigvec to use for each data type
if(update_c){
  c_single = map_dbl(Z0, function(z){
    estimateNumberOfClustersGivenGraph(z, NUMC = num_eig)[[1]]
  })
  print(c_single)
}


F0 = lapply(1:S, function(i){
  L = diag(rowSums(Z0[[i]]))-Z0[[i]]
  F_eig1 = eig1(L, c = c_single[i], isMax = 0)$eigvec
  # F_eig1 = dn.cimlr(F_eig1, "ave")
  F_eig1
})
Z_cur = Z0
F_cur = F0
w_cur = rep(1/S,S)
Y_cur = matrix(0, n, c)
# optimization ----

# Update Z ----
Z_pre = Z_cur
for(s in 1:S){
  # updata each data type separately
  F_eig1 = F_cur[[s]]
  distX = initial_list[[s]]$distX
  distX1 = initial_list[[s]]$distX1
  idx = initial_list[[s]]$idx
  lambda = initial_list[[s]]$lambda
  r = initial_list[[s]]$lambda
  distf = L2_distance_1(t(F_eig1),t(F_eig1))
  A = array(0,c(n,n))
  b = idx[,2:dim(idx)[2]]
  a = apply(array(0,c(n,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:n })
  inda = cbind(as.vector(a),as.vector(b)) # rank of each row aligned
  ad = (distX[inda]+lambda*distf[inda])/2/r
  dim(ad) = c(n,ncol(b))

  # call the c function for the optimization
  c_input = -t(ad)
  c_output = t(ad)
  ad = t(.Call("projsplx_R",c_input,c_output))
  A[inda] = as.vector(ad)
  A[is.nan(A)] = 0
  A = (A + t(A)) / 2
  # Z_cur[[s]] = (1 - eta) * Z_cur[[s]] + eta * A
  Z_cur[[s]] = A
  # if(network_diffusion){Z_cur[[s]] = network.diffusion(Z_cur[[s]], k)}
  # Z_cur[[s]] = dn(Z_cur[[s]],"ave")
}
# Update number of eigen-vector to use in single
if(update_c){
  c_single = map_dbl(Z_cur, function(z){
    estimateNumberOfClustersGivenGraph(z, NUMC = num_eig)[[1]]
  })
  print(c_single)
}

