lihc_surv_tib_nc6 = res_tib_integ_nc %>%
  filter(cancer == "lihc" & n_cluster == 6 & dist_type == "dist") %>%
  dplyr::select(cancer, kernel_type,method, surv_tib) %>%
  unnest(surv_tib)

surv_6 = ggsurvplot_facet(survfit(Surv(days_to_death, censoring)~cluster,
                                  data = lihc_surv_tib_nc6),
                          data = lihc_surv_tib_nc6,
                          facet.by = c("kernel_type", "method"))+
  ggtitle("Survival curve for lihc using integrated original distance, nc = 6")

surv_6

part_cim_kernel = lihc_surv_tib_nc6 %>%
  filter(method == "part_cimlr" & kernel_type == "kernel")
table(part_cim_kernel$cluster)
