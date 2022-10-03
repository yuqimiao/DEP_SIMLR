# Header --
setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:8000
n = N[task_ID]

library(abSNF)
library(tidyverse)
dist2_w_n = function (X, C, weight)
{
  for (i in 1:dim(X)[1]) {
    X[i, ] = sqrt(weight) * X[i, ]
    C[i, ] = sqrt(weight) * C[i, ]
  }
  ndata = nrow(X)
  ncentres = nrow(C)
  sumsqX = rowSums(X^2)
  sumsqC = rowSums(C^2)
  XC = 2 * (X %*% t(C))
  res = matrix(rep(sumsqX, times = ncentres), ndata, ncentres) +
    t(matrix(rep(sumsqC, times = ndata), ncentres, ndata)) -
    XC
  res[res < 0] = 0
  return(res)
}
# read in ----
dir = "data_20220308"
name_tib = expand_grid(cancer = c("brca","kirp","lihc"),
                       type = c("ge","me")) %>%
  mutate(id = seq_along(cancer))

cancer = name_tib[n,]$cancer
type = name_tib[n,]$type
cancer_path = paste(dir, cancer,"cancer", paste(cancer,type,"cancer.Rdata", sep = "_"), sep = "/")
normal_path = paste(dir, cancer,"normal", paste(cancer,type,"normal.Rdata", sep = "_"), sep = "/")
cancer_data_name = load(cancer_path)
normal_data_name = load(normal_path)
cancer_data = eval(parse(text = cancer_data_name))
normal_data = eval(parse(text = normal_data_name))

cancer_data = cancer_data[,colnames(cancer_data)!="NAME"]

# normalization
library(SNFtool)
cancer_data = standardNormalization(cancer_data)
normal_data = standardNormalization(normal_data)

## bind samples ----
all_data = cbind(normal_data, cancer_data)
group = c(rep(0,dim(normal_data)[2]),rep(1,dim(cancer_data)[2])) # 0 for normal, 1 for case

# dim(all_data)
# length(group)

## 2 sample t test ----

log10_p_value_vec = map_dbl(1:dim(all_data)[1], function(i){
  print(i)
  x = all_data[i,]
  y = group
  control = x[y==0]
  case = x[y==1]
  case = as.numeric(case[!is.na(case)])
  control = as.numeric(control[!is.na(control)])
  if(length(unique(case)) == 1 & length(unique(control)) == 1){
    if(unique(case) == unique(control)){
      p = 1
    }else{
      p = broom::tidy(t.test(case,control))$p.value
    }
  }else{
    p = broom::tidy(t.test(case,control))$p.value
  }
  # p = p + .Machine$double.eps
  log10(p)
})

log10_p_value_vec[is.infinite(log10_p_value_vec)] = min(log10_p_value_vec[!is.infinite(log10_p_value_vec)])
w_vec = -log10_p_value_vec/sum(-log10_p_value_vec)

library(abSNF)
dist = dist2(t(cancer_data),t(cancer_data))
dist_w = dist2_w_n(t(cancer_data),t(cancer_data),w_vec)

saveRDS(dist, file = paste(dir, "/output/", cancer_data_name, "_dist.rds", sep = ""))
saveRDS(dist_w, file = paste(dir, "/output/", cancer_data_name, "_dist_w.rds", sep = ""))
saveRDS(w_vec, file = paste(dir, "/output/", cancer_data_name, "_w.rds", sep = ""))

