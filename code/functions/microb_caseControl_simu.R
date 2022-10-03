# data refered ----
# source("code/local_header.R")
data(throat.tree); data(dd)
OTUnames =throat.tree$tip.label =  paste0('OTU', throat.tree$tip.label) ## change OTU names, also need to change in distance calc
names(dd$pi) = paste0('OTU', names(dd$pi)); names(dd$theta) = paste0('OTU', names(dd$theta))

# do clustering ----
ncluster = 20
clustering = pam(as.dist(cophenetic(throat.tree)), ncluster, diss = TRUE)$clustering
p_est = dd$pi # mle of prob
p_est = p_est[names(clustering)]
p_clus = sort(tapply(p_est, clustering, sum), decreasing = T) # abundance level of the 20 clusters

# some parameters for distribution of OTU counts ----
theta = dd$theta; gplus = (1-theta)/theta; g_est = p_est*gplus

# final sig_tib ----
# source("code/local_header.R")
# source("code/weekly_analysis_script/Week2021-10-07/related_OTUselection.R")
# source("code/weekly_analysis_script/Week2021-10-07/unrelated_OTUselection.R")
# source("code/weekly_analysis_script/Week2021-10-07/related_unrelated_OTUselection.R")
#
# sig_tib = rbind(sig_tib, mix_sig)
#
# save(sig_tib, file = "/Volumes/Yuqi_Backup/Projects/Research_data/multiomics_integretion/microbiome data/simu_otu_tib.Rdata")

# parameters ----

load("/Volumes/Yuqi_Backup/Projects/Research_data/multiomics_integretion/microbiome data/simu_otu_tib.Rdata")
microb_sim_2 = function(sig_tib,
                     eta_perc = 0.0015, # right tail
                     # eta_thresh = rep(3, 3)
                     eta_thresh = NA,
                     nsam = 100000,
                     b = 1,
                     nsam_perc = 100,
                     type_cur = 1,
                     plot = F){
  sim = sim_data_g_est(g_est = g_est, nSam = nsam)

  sim_tib = sig_tib %>%
    filter(type == type_cur) %>%
    nest(data = -type) %>%
    mutate(sample_tib = map(data, function(data){
      eta_ls = list()
      ls = map(data$info_otus,"otu")
      names(ls) = data$effect
      for(i in 1:length(ls)){
        B = rep(b, length(ls[[i]]))
        if(str_detect(names(ls)[[i]], "abun")){
          sim_cur = sim/rowSums(sim)
        }else{
          sim_cur = sim
          sim_cur[sim_cur>0]=1
        }
        X = sim_cur[,ls[[i]]]
        eta_cur = scale(rowSums(X%*%B))
        if(plot==T){hist(eta_cur, main = paste("eta", i, sep = " "))}
        eta_ls[[paste("eta",i,sep = "")]] = eta_cur

      }
      eta_tib = as_tibble(eta_ls) %>% mutate(id = seq_along(eta1))
      if(is.na(eta_thresh)){eta_thresh = apply(eta_tib[,1:3], 2, function(x) quantile(x, (1-eta_perc)))}
      print(eta_thresh)
      eta_tib = eta_tib %>%
        mutate(H = pmap(list(e1 = eta1, e2 = eta2, e3 = eta3),
                        function(e1,e2,e3){
                          tibble(h1 = e1>=eta_thresh[1],
                                 h2 = e2>=eta_thresh[2],
                                 h3 = e3>=eta_thresh[3])
                        })
        )
      eta_tib = eta_tib %>%
        unnest(H) %>%
        mutate(C = pmap_chr(list(h1 = h1, h2 = h2, h3 = h3),
                            function(h1, h2, h3){
                              paste(which(c(h1,h2,h3)), collapse = ",")
                            })) %>%
        # filter(C!="" & !str_detect(C, ",")) %>%
        filter(C!="") %>%
        mutate(C = ifelse(str_detect(C, ","), str_split(C, ",")[[1]][1], C)) %>%
        nest(data = -C) %>%
        mutate(samples = map2(data,C, function(data, c){
          dat = data[,c("id", paste("eta",c, sep = ""))]
          colnames(dat)[2] = "eta"
          dat = dat %>% arrange(eta)
          dat[1:nsam_perc, "id"]
        })) %>% dplyr::select(-data)
    }))

  sim_tib = sim_tib %>% dplyr::select(-data) %>% unnest(sample_tib) %>% unnest(samples)
  sim = sim[sim_tib$id,]
  dist_list = dist_ls(sim, OTUnames = colnames(sim))

  return(list(sim =sim,
              truelabel = as.numeric(sim_tib$C),
              dist_list = dist_list))

}
# sim data ----

# tmp = microb_sim_2(sig_tib = sig_tib, nsam = 1000, nsam_perc = 10)
