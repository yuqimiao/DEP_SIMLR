# dependencies
library(tidyverse)
library(SNFtool)

library(flexsurv)
library(survival)
library(MASS)
library(survminer)


cluster_tib_dir = "/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap"
clinical_dir = "/Users/miaoyuqi/Desktop/temp_data/clinical"
cluster_path_ls = list.files(cluster_tib_dir)
cluster_path_ls = cluster_path_ls[str_detect(cluster_path_ls, "final_clust_tib.rds")]
clinical_path_ls = list.files(clinical_dir)
recipe = tibble(cluster_tib_path = cluster_path_ls,
                cancer = str_sub(cluster_path_ls, 1, 4)) %>%
  mutate(cluster_tib_path = paste(cluster_tib_dir, "/",cluster_path_ls, sep = "")) %>%
  mutate(clinical_path = map_chr(cancer, function(cancer){
    print(cancer)
    mask = str_detect(clinical_path_ls, regex(cancer, ignore_case = T))
    print(clinical_path_ls[mask])
    paste(clinical_dir, "/", clinical_path_ls[mask], sep = "")
  }))

surv_tib_all = NULL
for(i in 1:nrow(recipe)){
  cur_cancer = recipe$cancer[[i]]
  clinical_path = recipe$clinical_path[[i]]
  cluster_tib_path = recipe$cluster_tib_path[[i]]

  print(cur_cancer)
  cluster_tib = readRDS(cluster_tib_path)
  common_sample = colnames(cluster_tib$sim_mat[[1]])
  clinical = read_delim(clinical_path, delim = ",")
  clinical = clinical[3:nrow(clinical),]
  survival_tib = clinical %>%
    mutate(id = str_sub(bcr_patient_barcode,-4,-1))%>%
    dplyr::select(id, days_to_death, days_to_last_follow_up) %>%
    mutate(days_to_death = as.numeric(days_to_death),
           days_to_last_follow_up = as.numeric(days_to_last_follow_up)) %>%
    mutate(censoring = map2_dbl(days_to_death, days_to_last_follow_up, function(dd, dlfu){
      if(is.na(as.numeric(dd))){
        censoring = ifelse(is.na(as.numeric(dlfu)),NA,0)
      }else{
        censoring = 1
      }
      return(censoring)
    })) %>%
    filter(!is.na(censoring)) %>%
    mutate(days_to_death = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death))

  surv_clust_tib = cluster_tib %>%
    dplyr::select(-c(data,sim_mat)) %>%
    nest(data = -c(cancer, data_type,method, kernel_type,c)) %>%
    mutate(tib = map(data, function(dat){
      tibble(id = str_sub(common_sample,1,4),
             cluster = dat$cluster[[1]]) %>%
        left_join(survival_tib, by = "id")
    })) %>%
    mutate(surv_pv = map(tib, function(tib){
      diff = survdiff(Surv(days_to_death, censoring)~cluster, data = tib)
      pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
    }))

  saveRDS(surv_clust_tib, file = paste(cluster_tib_dir,"/", cur_cancer,"_surv_clust_tib.rds", sep = "" ))
  surv_tib_all = rbind(surv_tib_all, surv_clust_tib)
}

saveRDS(surv_tib_all, file = paste(cluster_tib_dir, "/all_surv_clust_tib.rds", sep = "" ))






