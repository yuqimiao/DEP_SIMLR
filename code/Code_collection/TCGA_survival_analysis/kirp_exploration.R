library(igraph)
obj_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/cancer_object/obj_tib.rds")

# oull out the res_list for kirp ge
res_kirp_ge = (obj_tib %>% filter(cancer == "kirp" &type == "ge" & dist_type == "dist" ) %>% pull(object))[[1]]
names(res_kirp_ge)

# est_n_clust didn't change after diffusion
res_kirp_ge$est_nclust_kernel
res_kirp_ge$est_nclust_diff_kernel

# table for the cluster
table(res_kirp_ge$cluster_kernel)
table(res_kirp_ge$cluster_diff_kernel)
# cluster nmi, very low
compare(res_kirp_ge$cluster_kernel, res_kirp_ge$cluster_diff_kernel)

# check neighb of a single subject

neighbor_difference = function(kernel,diff_kernel,subj_id = 1){
  n = nrow(kernel)
  k = floor(sqrt(n))
  par(mfrow = c(2,1))
  plot(kernel[subj_id,], ylim = c(0.7,1))
  plot(diff_kernel[subj_id,], ylim = c(0.7,1))

  sort_kernel_vec = sort(kernel[subj_id,], index.return =T, decreasing = T)
  sort_diff_kernel_vec = sort(diff_kernel[subj_id,], index.return =T, decreasing = T)

  sort_kernel_vec$ix[1:k]
  sort_diff_kernel_vec$ix[1:k]


  par(mfrow = c(1,2))
  plot(x = 1:k, y = sort_kernel_vec$x[1:k],xaxt = "n", ylab = "kerbel", xlab = "index", ylim = c(0.7,1))
  axis(1, at = 1:k,labels = as.character(sort_kernel_vec$ix[1:k]),las = 2)
  abline(h=0.85, col = "lightgray")
  plot(x = 1:k, y = sort_diff_kernel_vec$x[1:k],xaxt = "n", ylab = "diff_kenrel", xlab = "index", ylim = c(0.7,1))
  axis(1, at = 1:k,labels = as.character(sort_diff_kernel_vec$ix[1:k]),las = 2)
  abline(h=0.85, col = "lightgray")

  return(list(sort_kernel_vec = sort_kernel_vec,
              sort_diff_kernel_vec = sort_diff_kernel_vec))
}

neighbor_difference(kernel = res_kirp_ge$kernel,
                    diff_kernel = res_kirp_ge$diff_kernel,
                    subj_id = 2)

