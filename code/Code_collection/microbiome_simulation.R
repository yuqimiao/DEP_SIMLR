
# 1 basic environment and packages [simulation_function.R] ----
library(tidyverse)
library('MiSPU'); library('GUniFrac')
## load data ----
data(throat.tree); data(dd)
OTUnames =throat.tree$tip.label =  paste0('OTU', throat.tree$tip.label) ## change OTU names, also need to change in distance calc
names(dd$pi) = paste0('OTU', names(dd$pi)); names(dd$theta) = paste0('OTU', names(dd$theta))

## do clustering ----
ncluster = 20
clustering = pam(as.dist(cophenetic(throat.tree)), ncluster, diss = TRUE)$clustering
p_est = dd$pi # mle of prob
p_est = p_est[names(clustering)]
p_clus = sort(tapply(p_est, clustering, sum), decreasing = T) # abundance level of the 20 clusters

## some parameters for distribution of OTU counts ----
theta = dd$theta; gplus = (1-theta)/theta; g_est = p_est*gplus

# 2 major function [microb_caseControl_simu.R] ----
load("/Volumes/Yuqi_Backup/Projects/Research_data/multiomics_integretion/microbiome data/simu_otu_tib.Rdata")
microb_sim_2 = function(sig_tib,
                        eta_perc = 0.0015, # right tail
                        # eta_thresh = rep(3, 3)
                        eta_thresh = NA,
                        nsam = 100000,
                        b = 1,
                        nsam_perc = 100,# number of samples per cluster
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
# 3 utils [simulation_function.R]----
set_info_OTUs = function(info_type = 3){
  ## load data ----
  data(throat.tree); data(dd)
  throat.tree$tip.label = paste0('OTU', throat.tree$tip.label)
  names(dd$pi) = paste0('OTU', names(dd$pi)); names(dd$theta) = paste0('OTU', names(dd$theta))

  ## do clustering ----
  ncluster = 20
  clustering = pam(as.dist(cophenetic(throat.tree)), ncluster, diss = TRUE)$clustering
  p_est = dd$pi # mle of prob
  p_est = p_est[names(clustering)]
  p_clus = sort(tapply(p_est, clustering, sum), decreasing = T) # abundance level of the 20 clusters

  ## some parameters for distribution of OTU counts ----
  theta = dd$theta; gplus = (1-theta)/theta; g_est = p_est*gplus
  ## Set info OTUs----
  if(type == 3){
    info_OTUs_ls = vector('list', length = 3) # list of info OTUs
    names(info_OTUs_ls) = c('Phy_related_abun', 'Phy_related_pres', 'Phy_unrelated_abun')

    # 1st set of info OTUs: the 3rd largest cluster. We will use abundance of these OTUs
    info_OTUs_ls[[1]] = names(which(clustering == names(p_clus)[3]))
    # length(info_OTUs_ls[[1]]) # 85 OTUs
    # sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1026

    # 2nd set of info OTUs: the 2nd largest cluster. We will use presence/absence of these OTUs
    info_OTUs_ls[[2]] = names(which(clustering == names(p_clus)[2]))
    # length(info_OTUs_ls[[2]]) # 57 OTUs
    # sum(p_est[info_OTUs_ls[[2]]]) # abundance level 0.1039

    # 3rd set of info OTUs: some phylo-unrelated OTUs with large abundances not in the previous 2 sets
    # We will use abundance of these OTUs
    temp_OTUs = names(sort(p_est, decreasing = T)[15:50])
    temp_TF = rep(T, length(temp_OTUs))
    for (j in 2:length(temp_OTUs)) { # to make sure OTUs are in different clusters
      if (clustering[temp_OTUs[j]] %in% clustering[temp_OTUs[1:(j-1)]]) { temp_TF[j] = F }
    }
    temp_OTUs = temp_OTUs[temp_TF]
    info_OTUs_ls[[3]] = temp_OTUs[-which(temp_OTUs %in% unlist(info_OTUs_ls[1:2]))][1:10]
    # length(info_OTUs_ls[[3]]) # 10 OTUs
    # sum(p_est[info_OTUs_ls[[3]]]) # abundance level 0.1006
    # length(unique(clustering[info_OTUs_ls[[3]]]))
    # sort(p_est[info_OTUs_ls[[3]]], decreasing = T)
  }else if(type == 1){
    info_OTUs_ls_Phy_related_abun = vector('list', length = 3) # list of info OTUs
    names(info_OTUs_ls) = c('Phy_related_abun_1', 'Phy_related_abun_2', 'Phy_related_abun_3')

    # 1st set of info OTUs: the 3rd largest cluster. We will use abundance of these OTUs
    info_OTUs_ls_Phy_related_abun[[1]] = names(which(clustering == names(p_clus)[2]))
    info_OTUs_ls_Phy_related_abun[[2]] = names(which(clustering == names(p_clus)[3]))
    info_OTUs_ls_Phy_related_abun[[3]] = names(which(clustering == names(p_clus)[4]))
    # length(info_OTUs_ls[[1]]) # 85 OTUs
    # sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1026
  }else if(type == 2){
    info_OTUs_ls = vector('list', length = 3) # list of info OTUs
    names(info_OTUs_ls) = c('Phy_related_pre/abs_1', 'Phy_related_pre/abs_2', 'Phy_related_pre/abs_3')

    # 1st set of info OTUs: the 3rd largest cluster. We will use abundance of these OTUs
    info_OTUs_ls_Phy_related_abun[[1]] = names(which(clustering == names(p_clus)[18]))
    info_OTUs_ls_Phy_related_abun[[2]] = names(which(clustering == names(p_clus)[19]))
    info_OTUs_ls_Phy_related_abun[[3]] = names(which(clustering == names(p_clus)[20]))
    # length(info_OTUs_ls[[1]]) # 85 OTUs
    # sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1026
  }

  return(info_OTUs_ls)
}
sim_data_g_est = function(nSam = 100, mu = 1000, size = 25, g_est)
{
  comm = matrix(0, nSam, length(g_est))
  comm.p = comm
  rownames(comm) = 1:nrow(comm)
  colnames(comm) = names(g_est)
  nSeq = rnbinom(nSam, mu = mu, size = size)
  for (i in 1:nSam) {
    comm.p[i, ] = MiSPU::rdirichlet(1, g_est)[1, ]
    comm[i, ] = rmultinom(1, nSeq[i], prob = comm.p[i, ])[,1]
  }
  return(comm)
}
UniFrac_BC = function(OTUtab = OTUtab) {
  OTU_freq_tab = t(apply(OTUtab, 1, function(x) x/sum(x)))
  return(as.matrix(stats::dist(OTU_freq_tab, method = 'manhattan'))/2)
}
UniFrac_HM = function(OTUtab = OTUtab) {
  OTUtab[OTUtab>0] = 1
  return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/ncol(OTUtab))
}
dist_ls = function(dataset,OTUnames, subtree = throat.tree){
  if(is.na(OTUnames)){OTUnames = colnames(dataset)}
  OTUtab = dataset[, OTUnames]
  unifracs = GUniFrac::GUniFrac(otu.tab = OTUtab, tree = subtree, alpha = 1)$unifracs
  distmat_ls = list(UniFrac_BC(OTUtab), UniFrac_HM(OTUtab), unifracs[,,'d_1'], unifracs[,,'d_UW'])
  names(distmat_ls) = c('BC', 'HM', 'Weighted', 'Unweighted')
  distmat_ls
}

# simulation ----
tmp = microb_sim_2(sig_tib = sig_tib,
                   type = 2,
                   eta_perc = 0.5,
                   nsam = 10000,
                   nsam_perc = 50)

# 0.0015 -- many 0
# 0.1 -- still many 0 but 1 is siginificantly higher

