# dependencies
library(tidyverse)
library(SNFtool)

library(flexsurv)
library(survival)
library(MASS)
library(survminer)
# functions
# The goal of the function is to tidy the clinical data
# input:
# 1. cur cancer: the current cancer we focused on, not used in the function process but output as a tibble indicator
# 2. clinical path: path where the cinilical data of the cancer are stored, need to have, download by:
# library(TCGAbiolinks)
# kirp_clinical=GDCquery_clinic("TCGA-KIRP", type = "clinical", save.csv = T)
# 3. cluster_res: cluster result of patients, need to have the names of the cluster res to be the last 4 character in the TCGA bar code

# output:
# survival_tib with sample index, cluster, days_to_death, days_to_last_follow_up, censoring
survival_tib_tidy = function(cur_cancer = "kirp",
                             clinical_path = "/Users/miaoyuqi/Desktop/temp_data/TCGA-KIRP_clinical.csv",
                             cluster_res){
  # extract survival clinical data ----
  clinical = read_delim(clinical_path, delim = ",")
  clinical_features = names(clinical)
  clinical_features_days = clinical_features[str_detect(clinical_features,"days")]

  survival_tib = clinical %>%
    dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_follow_up) %>%
    mutate(censoring = map2_dbl(days_to_death, days_to_last_follow_up, function(dd, dlfu){
      if(is.na(dd)){
        censoring = ifelse(is.na(dlfu),NA,0)
      }else{
        censoring = 1
      }
      return(censoring)
    })) %>%
    filter(!is.na(censoring)) %>%
    mutate(days_to_death = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)) %>%
    mutate(sample_ind = str_sub(bcr_patient_barcode, -4, -1)) %>%
    dplyr::select(sample_ind, days_to_death, days_to_last_follow_up, censoring)

  survival_tib = left_join(tibble(sample_ind = names(cluster_res),
                                  cluster = cluster_res),
                           survival_tib,
                           by = "sample_ind")
   return(survival_tib)
  }

# data read in
obj_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/cancer_object/obj_tib.rds")
kernel_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/kernel_tib.rds")
benchmark_tib = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/benchmark_tib.rds")
clinical_dir = "/Users/miaoyuqi/Desktop/temp_data"

# result tidy ----

# add survival to integrated result tib, with estimated number of clusters of summed similarity  ----
res_tib_integ = benchmark_tib %>% # tidy, make the columns of res to a long column
  pivot_longer(c(SNF,part_cimlr, cimlr),
               names_to = "method",
               values_to = "res") %>%
  mutate(sample_name = map(kernel_list, function(ls) colnames(ls[[1]]))) %>%
  mutate(cluster_similarity = map2(res, sample_name, function(res, name){
    S = res$S
    cluster = res$cluster
    names(cluster) = name
    colnames(S) = rownames(S) = name
    tibble(cluster = list(cluster),
           similarity = list(S))
  })) %>% unnest(cluster_similarity) %>%
  dplyr::select(-c(kernel_list, res, est_nclust_single))

res_tib_integ = res_tib_integ %>% # addTCGA result
  mutate(surv_tib = map2(cancer, cluster, function(cancer, cluster){
    clinical_path = paste(clinical_dir, "/TCGA-", str_to_upper(cancer),"_clinical.csv", sep = "")
    print(clinical_path)
    survival_tib_tidy(cur_cancer = cancer,
                      clinical_path = clinical_path,
                      cluster_res = cluster)
  }))%>% # add survival _pv
  mutate(surv_pv = map(surv_tib, function(tib){
    diff = survdiff(Surv(days_to_death, censoring)~cluster, data = tib)
    pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
  }))
# 2 cancer * 2 dist_type * 6 methods (2 kenrel type * 3 methods, concatenated)
saveRDS(res_tib_integ, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_integ.rds")

# add survival to single data result tib ----
res_tib_single = obj_tib %>%
  mutate(cluster_similarity = map(object, function(obj){
    tibble(kernel_type = c("kernel", "diff_kernel"),
           est_nclust = c(obj$est_nclust_kernel[[1]],obj$est_nclust_diff_kernel[[1]]),
           cluster = list(obj$cluster_kernel, obj$cluster_diff_kernel),
           similarity = list(obj$kernel, obj$diff_kernel))
  })) %>% unnest(cluster_similarity) %>%
  dplyr::select(cancer, type, dist_type, est_nclust, kernel_type, cluster, similarity)

res_tib_single = res_tib_single %>%
  mutate(surv_tib = map2(cancer, cluster, function(cancer, cluster){
    clinical_path = paste(clinical_dir, "/TCGA-", str_to_upper(cancer),"_clinical.csv", sep = "")
    print(clinical_path)
    survival_tib_tidy(cur_cancer = cancer,
                      clinical_path = clinical_path,
                      cluster_res = cluster)
  })) %>% # add survival _pv
  mutate(surv_pv = map(surv_tib, function(tib){
    diff = survdiff(Surv(days_to_death, censoring)~cluster, data = tib)
    pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
  }))
saveRDS(res_tib_single, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_single.rds")




# survival infor extraction ----
surv_pv = rbind(res_tib_single %>% # single
                  mutate(method = "single") %>%
                  dplyr::select(-c(cluster, similarity, surv_tib)) %>%
                  unnest(surv_pv),

                res_tib_integ %>% # integration
                  mutate(type = "integration") %>%
                  dplyr::select(-c(cluster, similarity, surv_tib, sample_name)) %>%
                  unnest(surv_pv)) %>%
  dplyr::select(c(cancer, dist_type, type, kernel_type, method, est_nclust, surv_pv))


saveRDS(surv_pv, "code/Code_collection/TCGA_survival_analysis/surv_pv.rds")
saveRDS(surv_pv, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/surv_pv.rds")


# cluster distribution infor extraction ----

cluster_dist_tib_single = res_tib_single %>%
  mutate(cluster_dist = map(cluster, function(c){
    tib = as_tibble(table(c))
    names(tib)[1] ="label"
    tib})) %>%
  unnest(cluster_dist) %>%
  mutate(method = "single") %>%
  dplyr::select(cancer, type, dist_type, kernel_type, method, est_nclust, label, n)


cluster_dist_tib_integ = res_tib_integ %>%
  mutate(cluster_dist = map(cluster, function(c){
    tib = as_tibble(table(c))
    names(tib)[1] ="label"
    tib})
    ) %>%
  unnest(cluster_dist) %>%
  mutate(type = "integration")%>%
  dplyr::select(cancer, type, dist_type, kernel_type, method, est_nclust, label, n)

cluster_dist_tib = rbind(cluster_dist_tib_single, cluster_dist_tib_integ)
saveRDS(cluster_dist_tib, "code/Code_collection/TCGA_survival_analysis/cluster_distribution.rds")
saveRDS(cluster_dist_tib, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/cluster_distribution.rds")

# add survival to integrated result tib, with tunning number of clusters ----
benchmark_nc_tunning = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/benchmark_nc_tunning.rds")
res_tib_integ_nc = benchmark_nc_tunning %>% # tidy, make the columns of res to a long column
  pivot_longer(c(SNF,part_cimlr, cimlr),
               names_to = "method",
               values_to = "res") %>%
  mutate(sample_name = map(kernel_list, function(ls) colnames(ls[[1]]))) %>%
  mutate(cluster_similarity = map2(res, sample_name, function(res, name){
    S = res$S
    cluster = res$cluster
    names(cluster) = name
    colnames(S) = rownames(S) = name
    tibble(cluster = list(cluster),
           similarity = list(S))
  })) %>% unnest(cluster_similarity) %>%
  dplyr::select(-c(kernel_list, res, est_nclust_single))

res_tib_integ_nc = res_tib_integ_nc %>% # addTCGA result
  mutate(surv_tib = map2(cancer, cluster, function(cancer, cluster){
    clinical_path = paste(clinical_dir, "/TCGA-", str_to_upper(cancer),"_clinical.csv", sep = "")
    print(clinical_path)
    survival_tib_tidy(cur_cancer = cancer,
                      clinical_path = clinical_path,
                      cluster_res = cluster)
  }))%>% # add survival _pv
  mutate(surv_pv = map(surv_tib, function(tib){
    diff = survdiff(Surv(days_to_death, censoring)~cluster, data = tib)
    pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
  }))

saveRDS(res_tib_integ_nc, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_integ_nc.rds")

surv_pv_integ_nc = res_tib_integ_nc%>%
  mutate(type = "integration") %>%
  dplyr::select(c(cancer, dist_type, type, kernel_type, method, n_cluster, surv_pv))
saveRDS(surv_pv_integ_nc, "/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/surv_pv_integ_nc.rds")

 # part_cimlr, kernel, c = 6 get 9e-10?


















