library(tidyverse)
source("code/functions/visualization_functions.R")
# tidy the input data and rand by weight, extract the top n_feat
# input
# feature_path: stored cancer feature data, row as feature, columns as sample, with column/row names
# common_sample_path, stored common samples with gene expression, mutation and methylation data, list with kirp, brca, lihc
# feature_weight_path: weight generated for the features based on case/control, no mutation data
# n_feat

# output:
# data: n_feat*n_common_sample_cur
# feature_weight_sort: sorted first 50 feature weights
top_feat_data_tidy = function(feature_path = "/Volumes/sw2206/yuqi/data_20220308/kirp/cancer/kirp_ge_cancer.Rdata",
                              common_sample_path = "/Volumes/sw2206/yuqi/data_20220308/common_sample.rds",
                              feature_weight_path = "/Volumes/sw2206/yuqi/data_20220308/output/kirp_ge_cancer_w.rds",
                              n_feat = 50){
  common_sample = readRDS(common_sample_path)
  data_name = load(feature_path)
  data = eval(parse(text = data_name))
  cancer = str_sub(data_name,1,4)
  # intersect with common sample
  common_sample_cur = common_sample[[which(str_sub(names(common_sample),1,4) == cancer)]]
  colnames(data) = str_sub(colnames(data), -4,-1)
  data = data[,common_sample_cur]

  # get weight ranked top 50
  fw = readRDS(feature_weight_path)
  names(fw) = rownames(data)
  fw_sort = sort(fw, decreasing = T)[1:n_feat]
  feat_sorted_names = names(fw_sort)
  data = data[feat_sorted_names,]

  return(list(data = data,
              feat_weight_sort = fw_sort))
}




# result read in
res_tib_integ_nc = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/res_tib_integ_nc.rds")
surv_pv_integ_nc = readRDS("/Volumes/sw2206/yuqi/data_20220308/output/benchmark_res/surv_pv_integ_nc.rds")

kirp_nc4 = res_tib_integ_nc %>%
  filter(cancer == "kirp" & n_cluster == 4) %>%
  filter(kernel_type == "diff_kernel" & dist_type == "dist")

surv_pv_integ_nc %>% filter(cancer == "kirp" & n_cluster == 4) %>%
  filter(kernel_type == "diff_kernel" & dist_type == "dist") %>%
  unnest(surv_pv)

kirp_surv_nc4 = res_tib_integ_nc %>%
  filter(dist_type == "dist" & cancer == "kirp"& n_cluster == 4) %>%
  dplyr::select(cancer, method, kernel_type,surv_tib) %>%
  unnest(surv_tib)

ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                         data = kirp_surv_nc4),
                 data = kirp_surv_nc4,
                 facet.by = c("kernel_type", "method"),
                 pval = T)+
  ggtitle("Survival curve for kirp using original distance, nc = 4")

# color bar ----

cluster_ls = kirp_nc4$cluster
names(cluster_ls) = kirp_nc4$method
g_color_bar = color_bar_cl(cluster_ls = cluster_ls)

# feature heatmap ----
sample_order = (tibble(sample_id = names(cluster_ls[[1]]),
                      cluster = cluster_ls$part_cimlr) %>%
  arrange(cluster) %>%
  pull(sample_id))

kirp_ge_ls = top_feat_data_tidy(feature_path = "/Volumes/sw2206/yuqi/data_20220308/kirp/cancer/kirp_ge_cancer.Rdata",
                                common_sample_path = "/Volumes/sw2206/yuqi/data_20220308/common_sample.rds",
                                feature_weight_path = "/Volumes/sw2206/yuqi/data_20220308/output/kirp_ge_cancer_w.rds",
                                n_feat = 100)
mat_ge = kirp_ge_ls$data
g_ge = feature_heatmap(mat = mat_ge,
                column_id_order = sample_order,
                row_id_order = rownames(data),
                name = "KIRP Expression top 100 heatmap",
                discrete = F)

kirp_me_ls = top_feat_data_tidy(feature_path = "/Volumes/sw2206/yuqi/data_20220308/kirp/cancer/kirp_me_cancer.Rdata",
                                common_sample_path = "/Volumes/sw2206/yuqi/data_20220308/common_sample.rds",
                                feature_weight_path = "/Volumes/sw2206/yuqi/data_20220308/output/kirp_me_cancer_w.rds",
                                n_feat = 100)
mat_me = kirp_me_ls$data
g_me = feature_heatmap(mat = mat_me,
                       column_id_order = sample_order,
                       row_id_order = rownames(mat_me),
                       name = "KIRP Methylation top 100 heatmap",
                       discrete = F)

ggarrange(plotlist = list(g_ge,g_color_bar),
          nrow = 2,
          label.x = "sample",
          align = "v",
          heights = c(4,2))
ggsave("code/Code_collection/TCGA_survival_analysis/kirp_gene_heatmap.png")
ggarrange(plotlist = list(g_me,g_color_bar),
          nrow = 2,
          label.x = "sample",
          align = "v",
          heights = c(4,2))
ggsave("code/Code_collection/TCGA_survival_analysis/kirp_methylation_heatmap.png")
