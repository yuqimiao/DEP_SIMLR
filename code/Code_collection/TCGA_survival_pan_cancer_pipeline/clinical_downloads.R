library(TCGAbiolinks)
imputed_ls_dir = "/Volumes/sw2206/yuqi/data_20220308/imputed_list"
recipe = tibble(path = list.files(imputed_ls_dir)) %>%
  mutate(cancer = str_sub(path, 1,4)) %>%
  mutate(path = paste(imputed_ls_dir, path, sep = "/"))

for(i in 1:nrow(recipe)){
  cur = paste("TCGA-", recipe$cancer[[i]],sep = "")
  clinical=GDCquery_clinic(cur, type = "clinical", save.csv = T)
}
