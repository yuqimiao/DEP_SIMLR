library(tidyverse)
library(impute)
mut_dir = "/Volumes/Yuqi_Backup/HPC Data/mutation"
mutation_tib = tibble(mutation = list.files(mut_dir)[str_detect(list.files(mut_dir), ".Rdata")])
mutation_tib = mutation_tib %>%
  mutate(cancer = str_sub(mutation, 1,4),
         mutation = paste(mut_dir, mutation, sep = "/"))

methy_dir = "/Volumes/Yuqi_Backup/HPC Data/methylation"
methylation_tib = tibble(methylation = list.files(methy_dir)[str_detect(list.files(methy_dir), ".Rdata")])
methylation_tib = methylation_tib %>%
  mutate(cancer = str_sub(methylation, 1,4),
         methylation = paste(methy_dir, methylation, sep = "/"))

expr_dir = "/Volumes/Yuqi_Backup/HPC Data/expression"
expression_tib = tibble(expression = list.files(expr_dir))
expression_tib = expression_tib %>%
  mutate(cancer = str_sub(expression, 1,4),
         expression = paste(expr_dir, expression, sep = "/"))

pan_cancer_tib = mutation_tib %>%
  filter(!(cancer %in% c("HNSC","TGCT","UCEC")))%>%
  inner_join(methylation_tib, by = "cancer") %>%
  inner_join(expression_tib, by = "cancer") %>%
  dplyr::select(cancer, expression, methylation, mutation)

common_sample_list = list()
data_list_all = list()
imputed_dir = "/Volumes/Yuqi_Backup/HPC Data/imputed_list"
for(i in 1:nrow(pan_cancer_tib)){
  cancer = pan_cancer_tib$cancer[[i]]
  ex = pan_cancer_tib$expression[[i]]
  me = pan_cancer_tib$methylation[[i]]
  mu = pan_cancer_tib$mutation[[i]]

  print(cancer)
  expr = read.table(ex, header = F)
  expr_sample = as.character(expr[[1]])
  expr_sample_nafilter = expr_sample[which(apply(expr,1,function(x) sum(is.na(x))/length(x) <=0.3))]

  expr_t = apply(expr[,-1],1, as.numeric)
  colnames(expr_t) = expr_sample


  load(me)
  me_sample = colnames(tumor)
  me_sample_nafilter = me_sample[which(apply(tumor,2,function(x) sum(is.na(x))/length(x) <=0.3))]

  load(mu)
  mu_sample = rownames(driver)
  mu_sample_nafilter = mu_sample[which(apply(driver,1,function(x) sum(is.na(x))/length(x) <=0.3))]

  sample_list = list(expr_sample, me_sample, mu_sample)
  common_sample = intersect(intersect(str_sub(expr_sample_nafilter,4,15), str_sub(me_sample_nafilter,4,15)), str_sub(mu_sample_nafilter,4,15))
  data_list = list(expr = expr_t, me = tumor, mu = t(driver)) # row as features, columns as samples

  data_list = lapply(data_list, function(dat){
    # subsect only common samples
    colnames(dat) = str_sub(colnames(dat), 4, 15)
    dat = dat[,common_sample]
    n = length(common_sample)
    # remove features with >30% missings
    dat = dat[which(apply(dat,1,function(x) sum(is.na(x))/length(x) <=0.3)),]

    # KNN imputation
    dat_imp = impute.knn(data = as.matrix(dat), k=sqrt(n))$data

    return(dat_imp)
  })
  saveRDS(data_list, file = paste(imputed_dir, "/", cancer,"_imputDat_ls.rds", sep = ""))
  common_sample_list[[cancer]] = common_sample
  data_list_all[[cancer]] = data_list
}

pan_cancer_tib = pan_cancer_tib %>%
  filter(!(cancer %in% c("HNSC","TGCT"))) %>%
  mutate(common_sample = common_sample_list) %>%
  mutate(n_sample = map_dbl(common_sample, function(x) length(x)))

saveRDS(pan_cancer_tib, "/Volumes/Yuqi_Backup/HPC Data/cancer_recipe.rds")

