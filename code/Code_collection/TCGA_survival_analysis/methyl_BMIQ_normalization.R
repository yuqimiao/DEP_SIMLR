setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
n = N[task_ID]

library(tidyverse)
library(minfi)
library(sva)
library(readr)
library(dplyr)
library(wateRmelon)
library(impute)

# file dir
# imputed_ls_dir = "/Volumes/Yuqi_Backup/HPC Data/imputed_list"
imputed_ls_dir = "data_20220308/imputed_list"

# imputed_ls recipe
recipe = tibble(path = list.files(imputed_ls_dir)) %>%
  mutate(cancer = str_sub(path, 1,4)) %>%
  mutate(path = paste(imputed_ls_dir, path, sep = "/"))

# methylation
# methy_dir = "/Volumes/Yuqi_Backup/HPC Data/methylation"
methy_dir = "data_20220308/methylation"
methylation_tib = tibble(methylation = list.files(methy_dir)[str_detect(list.files(methy_dir), ".Rdata")]) # tib containing the methylation file paths for each cancer
methylation_tib = methylation_tib %>%
  mutate(cancer = str_sub(methylation, 1,4),
         methylation = paste(methy_dir, methylation, sep = "/"))
hg38_ann = read_delim(paste(methy_dir,"/HM450.hg38.manifest.gencode.v36.tsv", sep = ""), delim = "\t")

methy_path_cur = methylation_tib$methylation[[n]]
cancer_cur = methylation_tib$cancer[[n]]
impls_path_cur = recipe$path[str_detect(recipe$path, cancer_cur)]

load(methy_path_cur)
imputed_ls = readRDS(impls_path_cur)
common_sample = colnames(imputed_ls$me)
colnames(tumor) = str_sub(colnames(tumor), 4, 16)

tumor = tumor[,common_sample]
grset_obj = makeGenomicRatioSetFromMatrix(mat = tumor, rownames = hg38_ann$probeID) # make the methylation beta matrix as a GenomicRatioSet object

# remove CpG sites with SNPs overlapping
grset_obj = dropLociWithSnps(grset_obj, snps=c("SBE","CpG"), maf=0)
methy_ann = as_tibble(getAnnotation(grset_obj))
methy = getBeta(grset_obj)
n = ncol(methy)

print(n)

# remove features with more than 30% missing
methy = methy[which(apply(methy,1,function(x) sum(is.na(x))/length(x) <=0.3)),]
# knn imputation
methy = impute.knn(methy, k = floor(sqrt(n)))$data

# perform BMIQ
print("start BMIQ")
methy_ann_cur = methy_ann %>% filter(Name %in% rownames(methy))
probes_type = ifelse(methy_ann_cur$Type == "I",1,2)
methy_norm = NULL
for(i in 1:ncol(methy)){
  print(paste("start sample", as.character(i)))
  tmp = BMIQ(methy[,i],design.v = probes_type,nfit = 500)
  # cur_column = as.data.frame(tmp$nbeta)
  cur_column = tmp$nbeta
  methy_norm = cbind(methy_norm, cur_column)
}

rownames(methy_norm) = rownames(methy)
colnames(methy_norm) = colnames(methy)
saveRDS(methy_norm, file = paste(methy_dir, "/processed_methy/",cancer_cur, "_filtered_imputed_BMIQ_methy.rds", sep = ""))

imputed_ls$me = methy_norm
saveRDS(imputed_ls, file = impls_path_cur)




