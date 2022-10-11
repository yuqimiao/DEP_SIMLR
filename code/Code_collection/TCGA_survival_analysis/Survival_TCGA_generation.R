# dependencies
source("code/functions/Partition_CIMLR_2.0.R")
source("code/functions/CIMLR_kernel_input.R")
library(SNFtool)
library(tidyverse)
# functions ----
single_cancer_obj_generate = function(distance_path,
                                      sigma = 2,
                                      alpha = 0.8,
                                      diffusion_form = "L1",
                                      output_path,
                                      common_sample){

  distance = readRDS(distance_path)
  sample_name = str_sub(rownames(distance),-4,-1 )
  n = nrow(distance)
  k = floor(sqrt(n))
  kernel = kernel_calculation(distance, k = k, sigma = sigma)
  diff_kernel = diffusion_enhancement(kernel, alpha = alpha, k = k, diffusion_form = diffusion_form)
  colnames(kernel)=rownames(kernel) = sample_name
  colnames(diff_kernel)=rownames(diff_kernel) = sample_name

  # For single data, we also only focus on those subjects with all omics data for comparison
  common_sample = as.character(common_sample)
  kernel = kernel[common_sample, common_sample]
  diff_kernel = diff_kernel[common_sample, common_sample]


  est_nclust_kernel = estimateNumberOfClustersGivenGraph(kernel, 2:11)
  est_nclust_diff_kernel = estimateNumberOfClustersGivenGraph(diff_kernel, 2:11)

  cluster_kernel = spectralClustering(kernel, est_nclust_kernel[[1]])
  cluster_diff_kernel = spectralClustering(diff_kernel, est_nclust_diff_kernel[[1]])
  names(cluster_kernel) = common_sample
  names(cluster_diff_kernel) = common_sample

  res_list = list(distance = distance,
                  kernel = kernel,
                  diff_kernel = diff_kernel,
                  est_nclust_kernel = est_nclust_kernel,
                  est_nclust_diff_kernel = est_nclust_diff_kernel,
                  cluster_kernel = cluster_kernel,
                  cluster_diff_kernel = cluster_diff_kernel)
  saveRDS(res_list, file = output_path)
  return("finish")
}

# Single data ----
## single data type object generation ----
dir = "data_20220308/output/"
common_sample = readRDS("data_20220308/common_sample.rds")
path_tib = expand_grid(cancer = c("kirp","lihc"),
                       type = c("ge","me", "mu"),
                       dist_type = c("dist","dist_w")) %>%
  filter(!(dist_type == "dist_w" & type == "mu")) %>%
  mutate(tag = paste(cancer,type,"cancer",dist_type, sep = "_")) %>%
  mutate(distance_path = paste(dir, tag, ".rds",sep = "")) %>%
  mutate(output_path = paste(dir, "cancer_object/", tag, "_obj", ".rds", sep = "")) %>%
  mutate(single_obj = pmap(list(d = distance_path, o = output_path, cancer = cancer), function(d,o, cancer){
    cs = common_sample[[paste(cancer, "_common_sample", sep = "")]]
    single_cancer_obj_generate(distance_path = d, output_path = o, common_sample = cs)
  }))
saveRDS(path_tib,"data_20220308/output/cancer_object/path_tib.rds")

obj_tib = path_tib %>% mutate(object = map(output_path, function(o){
  readRDS(o)
}))
saveRDS(obj_tib,"data_20220308/output/cancer_object/obj_tib.rds")






# integrated data ----


kernel_tib = obj_tib %>%
  nest(data = -c(cancer, dist_type)) %>%
  mutate(tib = map(data, function(dat){ # extract kernel
    kernel_ls = map(dat$object,"kernel")
    diff_kernel_ls = map(dat$object,"diff_kernel")
    common_sample = Reduce("intersect",lapply(kernel_ls, function(x) colnames(x)))
    kernel_ls = map(kernel_ls, function(x) x[common_sample,common_sample])
    diff_kernel_ls = map(diff_kernel_ls, function(x) x[common_sample,common_sample])
    kernel_nc = map_dbl(dat$object,function(o)o$est_nclust_kernel[[1]])
    diff_kernel_nc = map_dbl(dat$object,function(o)o$est_nclust_diff_kernel[[1]])
    tib = tibble(kernel_type = c("kernel", "diff_kernel"),
                 kernel_list = list(kernel_ls, diff_kernel_ls),
                 est_nclust_single = list(kernel_nc,diff_kernel_nc)
    )
    return(tib)
  })) %>%
  dplyr::select(cancer, dist_type, tib) %>%
  unnest(tib)

saveRDS(kernel_tib, file = "data_20220308/output/benchmark_res/kernel_tib.rds")

# SNF
benchmark_tib = kernel_tib %>%
  mutate(est_nclust = map_dbl(kernel_list, function(ls){
    ave = Reduce("+", ls)/length(ls)
    c = estimateNumberOfClustersGivenGraph(ave, 2:11)[[1]]
    return(c)
  })) %>%
  mutate(SNF = map2(kernel_list,est_nclust, function(ls,c){
    n = nrow(ls[[1]])
    k = floor(sqrt(n))
    S_all = SNF(ls, k)
    cluster = spectralClustering(S_all,c)
    return(list(S = S_all,
                cluster = cluster))
  }))
# add part_cim
benchmark_tib = benchmark_tib%>%
  mutate(part_cimlr = pmap(list(ls = kernel_list, c_single = est_nclust_single, c = est_nclust), function(ls,c_single,c){
    n = nrow(ls[[1]])
    k = floor(sqrt(n))
    print("here")
    res = part_cimlr(kernel_list = ls,k = k, c = c,neig_single = c_single, update_neig = F)
    return(list(S = res$Y %*% t(res$Y),
                cluster = res$cluster))
  }))

# add cim
benchmark_tib = benchmark_tib%>%
  mutate(cimlr = pmap(list(ls = kernel_list, c_single = est_nclust_single, c = est_nclust), function(ls,c_single,c){
    n = nrow(ls[[1]])
    k = floor(sqrt(n))
    res = CIMLR_kernel(kernel_list = ls,  c = c, k = k)
    return(list(S = res$S,
                cluster = res$cluster))
  }))


saveRDS(benchmark_tib, "data_20220308/output/benchmark_res/benchmark_tib.rds")

# add cluster result with different number of clusters as options
kernel_tib = readRDS("data_20220308/output/benchmark_res/kernel_tib.rds")

benchmark_nc_tune_tib = tibble(n_cluster = 3:6) %>%
  mutate(res_tib = map(n_cluster, function(nc){
    print(nc)
    tib = kernel_tib %>%
      mutate(SNF = map(kernel_list, function(ls){
        cat("SNF", nc)
        n = nrow(ls[[1]])
        k = floor(sqrt(n))
        S_all = SNF(ls, k)
        cluster = spectralClustering(S_all, nc)
        return(list(S = S_all,
                    cluster = cluster))
      }))%>%
      mutate(part_cimlr = map2(kernel_list, est_nclust_single, function(ls,c_single){
        cat("part_cimlr", nc)
        n = nrow(ls[[1]])
        k = floor(sqrt(n))
        res = part_cimlr(kernel_list = ls, k = k, neig_single = c_single, c = nc, update_neig = F)
        return(list(S = res$Y %*% t(res$Y),
                    cluster = res$cluster))
      })) %>%
      mutate(cimlr = map(kernel_list, function(ls){
        n = nrow(ls[[1]])
        k = floor(sqrt(n))
        cat("cimlr", nc)
        res = CIMLR_kernel(kernel_list = ls,  c = nc, k = k)
        return(list(S = res$S,
                    cluster = res$cluster))
      }))
    saveRDS(tib, paste("data_20220308/output/benchmark_res/benchmark_nc_tunning_", nc, ".rds", sep = ""))
    return(tib)
  }))

benchmark_nc_tune_tib = benchmark_nc_tune_tib %>%
  unnest(res_tib) %>%
  dplyr::select(cancer, dist_type,kernel_type, kernel_list, n_cluster, est_nclust_single, SNF, part_cimlr, cimlr)

saveRDS(benchmark_nc_tune_tib, paste("data_20220308/output/benchmark_res/benchmark_nc_tunning.rds", sep = ""))



