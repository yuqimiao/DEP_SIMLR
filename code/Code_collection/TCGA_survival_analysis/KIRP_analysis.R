# On original distance
# original data generated from "/Users/miaoyuqi/研究/Shuang project/multiomics_integation/code/Code_collection/TCGA_survival_analysis/Survival_TCGA_analysis.R"
# Compare diffusion and partition on integration----

library(tidyverse)
library(SNFtool)

library(flexsurv)
library(survival)
library(MASS)
library(survminer)

clinical = read_delim("/Users/miaoyuqi/Desktop/temp_data/TCGA-KIRP_clinical.csv", delim = ",")
res_tib_integ = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_integ.rds")
res_tib_single = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_single.rds")
kirp_surv_tib = res_tib_integ %>%
  filter(dist_type == "dist") %>%
  dplyr::select(cancer, method, kernel_type,surv_tib) %>%
  unnest(surv_tib)

ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                         data = kirp_surv_tib),
                 data = kirp_surv_tib,
                 facet.by = c("kernel_type", "method"))+
  ggtitle("Survival curve for kirp using original distance")

ggsave("code/Code_collection/TCGA_survival_analysis/kirp_dist_surv.png",width = 10, height = 6)

kirp_surv_tib_ge = res_tib_single %>%
  filter(dist_type == "dist" & type == "ge") %>%
  dplyr::select(cancer, kernel_type,surv_tib) %>%
  unnest(surv_tib)

ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                         data = kirp_surv_tib_ge),
                 data = kirp_surv_tib_ge,
                 facet.by = c("kernel_type"))+
  ggtitle("Survival curve for kirp using original ge distance")


kirp_surv_tib_me = res_tib_single %>%
  filter(dist_type == "dist" & type == "me") %>%
  dplyr::select(cancer, kernel_type,surv_tib) %>%
  unnest(surv_tib)

ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                         data = kirp_surv_tib_me),
                 data = kirp_surv_tib_me,
                 facet.by = c("kernel_type"))+
  ggtitle("Survival curve for kirp using original me distance")


# effect of diffusion ----
kirp_surv_tib_w = res_tib_integ %>%
  filter(dist_type == "dist_w") %>%
  dplyr::select(cancer, method, kernel_type,surv_tib) %>%
  unnest(surv_tib)

ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                         data = kirp_surv_tib_w),
                 data = kirp_surv_tib_w,
                 facet.by = c("kernel_type", "method"))+
  ggtitle("Survival curve for kirp using weighted distance")

ggsave("code/Code_collection/TCGA_survival_analysis/kirp_dist_w_surv.png",width = 10, height = 6)


# effect from integration

kirp_surv_pv = surv_pv %>% filter(cancer == "kirp")
lihc_surv_pv = surv_pv %>% filter(cancer == "lihc")


# KIRP diffusion effect explore ----

library(igraph)
obj_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/cancer_object/obj_tib.rds")

# pull out the res_list for kirp ge
res_kirp_ge = (obj_tib %>% filter(cancer == "kirp" &type == "ge" & dist_type == "dist" ) %>% pull(object))[[1]]
names(res_kirp_ge)

# est_n_clust didn't change after diffusion
res_kirp_ge$est_nclust_kernel
res_kirp_ge$est_nclust_diff_kernel

# table for the cluster distribution
table(res_kirp_ge$cluster_kernel)
table(res_kirp_ge$cluster_diff_kernel)

# cluster nmi
compare(res_kirp_ge$cluster_kernel, res_kirp_ge$cluster_diff_kernel, "nmi")

# check neighb of a single subject

# input: 2 similarity mat, kernel, diff kernel and the subject id which we want to look at, the last 4 digits of the tcga bar code
# output: sorted kernel/diff_kernel vec for the subject, and the plot showing the neighbor difference
neighbor_difference = function(kernel,diff_kernel,subj_id){
  n = nrow(kernel)
  k = floor(sqrt(n))
  par(mfrow = c(2,1))
  plot(kernel[subj_id,], ylim = c(0.7,1))
  plot(diff_kernel[subj_id,], ylim = c(0.7,1))

  sort_kernel_vec = sort(kernel[as.character(subj_id),], index.return = T, decreasing = T)
  sort_diff_kernel_vec = sort(diff_kernel[as.character(subj_id),], index.return = T, decreasing = T)



  par(mfrow = c(1,2))
  plot(x = 1:k, y = sort_kernel_vec$x[1:k],xaxt = "n", ylab = "kernel", xlab = "index", ylim = c(0.7,1), main = paste("subject", subj_id, "neighbor change"))
  axis(1, at = 1:k,labels = as.character(colnames(kernel)[sort_kernel_vec$ix[1:k]]),las = 2)
  abline(h=0.85, col = "lightgray")
  plot(x = 1:k, y = sort_diff_kernel_vec$x[1:k],xaxt = "n", ylab = "diff_kernel", xlab = "index", ylim = c(0.7,1))
  axis(1, at = 1:k,labels = as.character(colnames(kernel)[sort_diff_kernel_vec$ix[1:k]]),las = 2)
  abline(h=0.85, col = "lightgray")

  return(list(sort_kernel_vec = sort_kernel_vec,
              sort_diff_kernel_vec = sort_diff_kernel_vec))
}

neighbor_difference(kernel = res_kirp_ge$kernel,
                    diff_kernel = res_kirp_ge$diff_kernel,
                    subj_id = "7286")

source("code/functions/Partition_CIMLR_2.0.R")
distance = res_kirp_ge$distance
kernel = res_kirp_ge$kernel
diff_kernel = res_kirp_ge$diff_kernel
cluster_kernel = res_kirp_ge$cluster_kernel
cluster_diff_kernel = res_kirp_ge$cluster_diff_kernel

# subset distance measure
colnames(distance) = rownames(distance) = str_sub(colnames(distance), -4,-1)
distance = distance[names(cluster_kernel), names(cluster_kernel)]

distance = dist_kernels(diff_kernel)

as_tibble(cmdscale(distance, 2))%>%
  mutate(kernel = cluster_kernel,
         diff_kernel = cluster_diff_kernel) %>%
  pivot_longer(cols = c(kernel, diff_kernel),
               names_to = "kernel_type",
               values_to = "cluster") %>%
  ggplot(aes(x = V1, y = V2, color = factor(cluster))) +
  geom_point()+
  facet_grid(.~kernel_type)

eig_kernel_GL = eigen(normalized_GL(kernel))
eig_diff_kernel_GL = eigen(normalized_GL(diff_kernel))

eig_kernel_GL$values[(nrow(kernel)-6):nrow(kernel)]
eig_diff_kernel_GL$values[(nrow(kernel)-6):nrow(kernel)]


# ID: 6846 survival detail ----
id_cur = "6846"
obj_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/cancer_object/obj_tib.rds")
# extract kernel and diff kernel for ge
# pull out the res_list for kirp ge
res_kirp_ge = (obj_tib %>% filter(cancer == "kirp" &type == "ge" & dist_type == "dist" ) %>% pull(object))[[1]]
res_kirp_me = (obj_tib %>% filter(cancer == "kirp" &type == "me" & dist_type == "dist" ) %>% pull(object))[[1]]


kernel_ge = res_kirp_ge$kernel
diff_kernel_ge = res_kirp_ge$diff_kernel

kernel_me = res_kirp_me$kernel
diff_kernel_me = res_kirp_me$diff_kernel

par(mfrow = c(2,1))
tmp1 = neighbor_difference(kernel_ge, diff_kernel_ge, subj_id = id_cur)
tmp2 = neighbor_difference(kernel_me, diff_kernel_me, subj_id = id_cur)
