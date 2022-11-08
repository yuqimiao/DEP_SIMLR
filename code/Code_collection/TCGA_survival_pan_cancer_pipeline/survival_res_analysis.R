library(tidyverse)
library(flexsurv)
library(survival)
library(MASS)
library(survminer)
source("code/functions/visualization_functions.R")
surv_tib_all = readRDS("/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/all_surv_clust_tib.rds")
## Cherry picking ----
surv_tib_check = surv_tib_all %>%
  mutate(sample_size = map_dbl(tib, function(tib) nrow(tib))) %>%
  mutate(min_clust_size = map_dbl(tib, function(tib){
    min(table(tib$cluster))
  })) %>%
  filter(min_clust_size > 10) %>%
  filter(data_type == "Integration") %>%
  dplyr::select(-c(data,tib)) %>%
  filter(surv_pv<=0.05) %>%
  unnest(surv_pv) %>%
  ungroup() %>%
  arrange(cancer,surv_pv)

kirc_tib = surv_tib_all %>% filter(cancer == "KIRC") %>%
  mutate(sample_size = map_dbl(tib, function(tib) nrow(tib))) %>%
  mutate(min_clust_size = map_dbl(tib, function(tib){
    min(table(tib$cluster))
  })) %>%
  # filter(min_clust_size > 20) %>%
  filter(c == 3 | method == "spectral") %>%
  unnest(surv_pv)

kirc_check = kirc_tib %>%
  dplyr::select(-c(data,tib))

luad_tib = surv_tib_all %>% filter(cancer == "LUAD") %>%
  mutate(sample_size = map_dbl(tib, function(tib) nrow(tib))) %>%
  mutate(min_clust_size = map_dbl(tib, function(tib){
    min(table(tib$cluster))
  })) %>%
  # filter(min_clust_size > 20) %>%
  filter(c == 4 | method == "spectral") %>%
  unnest(surv_pv)

luad_check = luad_tib %>%
  dplyr::select(-c(data,tib))

kirp_tib = surv_tib_all %>% filter(cancer == "KIRP") %>%
  mutate(sample_size = map_dbl(tib, function(tib) nrow(tib))) %>%
  mutate(min_clust_size = map_dbl(tib, function(tib){
    min(table(tib$cluster))
  })) %>%
  # filter(min_clust_size > 20) %>%
  filter(c == 4 | method == "spectral") %>%
  unnest(surv_pv)

kirp_check = kirp_tib %>%
  dplyr::select(-c(data,tib))


## visual_function: ----

# input:
# 1. surv_tib: cancer, data_type, method, kernel_type, tib(for cluster res)
# 2. imp_ls: imputed list with expr, me, mu matrix
# 3. feat_rank: list with feature ranking of expr, me, mu
# 4. n_feat: number of featuers to show in the heatmap
# Output:
# g_sing_surv: the single data type survival plot
# g_integ_surv: the integrated survival plot
# g_me: the me feature heatmap
# g_expr: the expr feature heatmap

cancer_visual = function(cancer,
                         surv_tib,
                         imp_ls,
                         feat_rank,
                         n_feat = 100){
  # survival plot ----
  g_integ_surv = ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                                          data = surv_tib),
                                  data = subset(surv_tib, data_type == "Integration"),
                                  facet.by = c("kernel_type", "method"))+
    ggtitle(paste("Integrated data survival curve for", cancer))

  g_sing_surv = ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                                         data = surv_tib),
                                 data = subset(surv_tib, data_type != "Integration"),
                                 facet.by = c("kernel_type", "data_type"))+
    ggtitle(paste("Single data survival curve for", cancer))

  # color bar ----
  clust_ls_tib = surv_tib %>%
    filter(data_type == "Integration" & kernel_type == "diff_kernel") %>%
    dplyr::select(-c(days_to_death, days_to_last_follow_up, censoring)) %>%
    nest(cluster = -c(cancer, data_type, method, kernel_type)) %>%
    mutate(cluster = map(cluster, function(tib){
      cluster = tib$cluster
      names(cluster) = tib$id
      cluster
    }))

  # Making cluster aligning plot
  clust_ls = clust_ls_tib$cluster
  names(clust_ls) = clust_ls_tib$method
  g_clust = color_bar_cl(clust_ls)+
    theme(legend.position = "top")

  # heatmap ----
  # change sample names to the 4 digits to match the cluster name
  imp_ls = lapply(imp_ls, name_matching)
  g_expr_heat = feature_heatmap(mat = imp_ls$expr,
                                column_id_order = names(sort(clust_ls$part_cimlr)),
                                row_id_order = feat_rank$expr$aggR[1:n_feat],
                                discrete = F,
                                xlab = "sample",
                                ylab = "feature",
                                name = NA)+
    theme(legend.position = "bottom")

  g_expr = ggarrange(plotlist = list(g_clust + ggtitle(paste(cancer, "expression")),
                                     g_expr_heat),
                     nrow = 2,
                     label.x = "sample",
                     align = "v",
                     heights = c(2,4))


  g_me_heat = feature_heatmap(mat = imp_ls$me,
                              column_id_order = names(sort(clust_ls$part_cimlr)),
                              row_id_order = rownames(imp_ls$me)[feat_rank$me$aggR[1:n_feat]],
                              discrete = F,
                              xlab = "sample",
                              ylab = "feature",
                              name = NA)+
    theme(legend.position = "bottom")

  g_me = ggarrange(plotlist = list(g_clust + ggtitle(paste(cancer, "methylation")),
                                   g_me_heat),
                   nrow = 2,
                   label.x = "sample",
                   align = "v",
                   heights = c(2,4))
  return(list(g_integ_surv = g_integ_surv,
              g_sing_surv = g_sing_surv,
              g_me = g_me,
              g_expr = g_expr
  ))
}

name_matching = function(x){
  colnames(x) = str_sub(colnames(x),1,4)
  if(is.null(rownames(x))){
    rownames(x) = 1:nrow(x)
  }
  x
}


# KIRP ----

# kirp_clust = readRDS("/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/KIRP_final_clust_tib.rds")
kirp_fr = readRDS("/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/KIRP_feature_rank.rds")
kirp_il = readRDS("/Volumes/sw2206/yuqi/data_20220308/imputed_list/KIRP_imputDat_ls.rds")
kirp_surv_tib = kirp_tib %>%
  dplyr::select(cancer, data_type, method, kernel_type, tib) %>%
  unnest(tib)

g_kirp = cancer_visual(cancer = "KIRP",
              surv_tib = kirp_surv_tib,
              imp_ls = kirp_il,
              feat_rank = kirp_fr)


saveRDS(g_kirp, file = "/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/kirp_plots.rds")


# KIRC ----

kirc_fr = readRDS("/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/KIRC_feature_rank.rds")
kirc_il = readRDS("/Volumes/sw2206/yuqi/data_20220308/imputed_list/KIRC_imputDat_ls.rds")
kirc_surv_tib = kirc_tib %>%
  dplyr::select(cancer, data_type, method, kernel_type, tib) %>%
  unnest(tib)

g_kirc = cancer_visual(cancer = "KIRC",
                       surv_tib = kirc_surv_tib,
                       imp_ls = kirc_il,
                       feat_rank = kirc_fr)
saveRDS(g_kirc, file = "/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/kirc_plots.rds")

# LUAD ----
luad_fr = readRDS("/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/LUAD_feature_rank.rds")
luad_il = readRDS("/Volumes/sw2206/yuqi/data_20220308/imputed_list/LUAD_imputDat_ls.rds")
luad_surv_tib = luad_tib %>%
  dplyr::select(cancer, data_type, method, kernel_type, tib) %>%
  unnest(tib)

g_luad = cancer_visual(cancer = "LUAD",
                       surv_tib = luad_surv_tib,
                       imp_ls = luad_il,
                       feat_rank = luad_fr)

saveRDS(g_luad, file = "/Volumes/sw2206/yuqi/data_20220308/pancancer_output/alpha_0.8_eigengap/picking_3_cancer/luad_plots.rds")
