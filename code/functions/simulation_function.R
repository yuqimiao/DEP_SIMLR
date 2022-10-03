# 1 multiomics(centroid-distributed)
## 1.1 function intro ----
# function: simulation_1
# perform 1 simulation
# input:
#   1. x -- sample size
#   2. eff_ratio -- % of effective variables within all features
#   3. eff_size -- the effective size of the effective variables
#   4. sub_ratio -- a vector of the ratio of each subtype within sample
# output:
#   GE and MI, rows are samples, columns are features

# variance mode is only for dive4: 1/2/3/4
# if variance mode is same, this means that the clusters are having same within-cluster variance
# 3 data type simulation
# adjust number of subtype based on given ratio

# variance_ratio is only for dive4: 1/2/3/4
# if variance_ratio is rep(1,4), this means that the clusters are having same within-cluster variance
# else if variance_ratio is 4 different numbers, this means the variance of each cluster is different with this ratio

library(stringr)
library(rlist)

## 1.2main function ----
simulation_3 = function(size = 200,sub_ratio =rep(0.25,4),eff_size = rep(5,3),sigma = rep(100,3), variance_ratio = list(data1 = rep(1,4), data2 = rep(1,4), data3 = rep(1,4)),
                        data_divide =c("12/34","13/24","14/23"),dist = rep("normal",3),
                        feature_num = rep(1000,3), signal_ratio = rep(0.05,3),
                        n1_ratio = 0.6,n2_ratio = 0.4,r1 = 1,r2 = 2,uninfo_r1 = 0,uninfo_r2 = 1){
  if(inherits(eff_size,"list")& length(eff_size) == 1){
    eff_size = rep(eff_size,3)
  }

  ## calculate index_cut for given subtype ratio
  index_cut = cumsum(sub_ratio)*size
  index = list()
  for(i in 1:length(index_cut)){
    if(i == 1){
      index[[i]] = 1:index_cut[1]
    }else{
      index[[i]] = (index_cut[[i-1]]+1):index_cut[[i]]
    }
  }

  data1 = get_data(size, index, data_divide[1], dist[1],eff_size[[1]],sigma[1],signal_ratio = signal_ratio[1],feature_num = feature_num[1],variance_ratio = variance_ratio$data1)
  data2 = get_data(size, index, data_divide[2], dist[2],eff_size[[2]],sigma[2],signal_ratio = signal_ratio[2],feature_num = feature_num[2],variance_ratio = variance_ratio$data2)
  data3 = get_data(size, index, data_divide[3], dist[3],eff_size[[3]],sigma[3],signal_ratio = signal_ratio[3],feature_num = feature_num[3],variance_ratio = variance_ratio$data3)

  weight1 = get_weight(n1_ratio = n1_ratio, n2_ratio = n2_ratio,
                       r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2,feature_num = feature_num[1],signal_ratio = signal_ratio[1])
  weight2 = get_weight(n1_ratio = n1_ratio, n2_ratio = n2_ratio,
                       r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2,feature_num = feature_num[2],signal_ratio = signal_ratio[2])
  weight3 = get_weight(n1_ratio = n1_ratio, n2_ratio = n2_ratio,
                       r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2,feature_num = feature_num[3],signal_ratio = signal_ratio[3])
  truelabel = NULL
  for (i in 1:length(index)){
    truelabel = c(truelabel, rep(i,length(index[[i]])))
  }
  return(list(data1 = data1, data2 = data2, data3 = data3, weight1 = weight1, weight2 = weight2, weight3 = weight3,truelabel = truelabel) )
}



simulation_2 = function(size = 150,sub_ratio =  rep(1/3,3),eff_size = rep(1,2),sigma = rep(2^2,2),data_divide =c("1/23","12/3"),dist = rep("normal",2),n1 = 50*0.6,n2 = 950*0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ## calculate index_cut for given subtype ratio
  index_cut = cumsum(sub_ratio)*size
  index = list()
  for(i in 1:length(index_cut)){
    if(i == 1){
      index[[i]] = 1:index_cut[1]
    }else{
      index[[i]] = (index_cut[[i-1]]+1):index_cut[[i]]
    }
  }

  data1 = get_data(size, index, data_divide[1], dist[1],eff_size[1],sigma[1])
  data2 = get_data(size, index, data_divide[2], dist[2],eff_size[2],sigma[2])

  weight1 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)
  weight2 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)

  truelabel = NULL
  for (i in 1:length(index)){
    truelabel = c(truelabel, rep(i,length(index[[i]])))
  }
  return(list(data1 = data1, data2 = data2, weight1 = weight1, weight2 = weight2,truelabel = truelabel) )
}


## 1.3 sub functions ----
get_data = function(size,
                    index,
                    data_divide,
                    dist,
                    eff_size,
                    sigma,
                    signal_ratio,
                    feature_num,
                    variance_ratio){
  data1 = matrix(0,size, feature_num)
  signal_num = feature_num*signal_ratio
  divide = lapply(X = str_split(str_split(data_divide,pattern = "/")[[1]], pattern = ""),
                  FUN = as.numeric
  )
  ind_list = list()
  for(i in 1: length(divide)){
    ind = NULL
    for(j in 1:length(divide[[i]])){
      ind = c(ind,index[[divide[[i]][[j]]]])
    }
    ind_list[[i]] = ind
  }

  if(length(eff_size) == 1 & length(ind_list) == 2){
      eff_size = c(-eff_size,eff_size)
  }

  if(length(eff_size) == 1 & length(ind_list) == 4){
    eff_size = c(-2*eff_size,-eff_size, eff_size, 2*eff_size)
  }

  if (dist == "normal"){
    for(i in 1:length(ind_list)){
      data1[ind_list[[i]], 1:signal_num] = rnorm(length(ind_list[[i]])*signal_num,
                                        eff_size[i],
                                        sqrt(sigma*variance_ratio[i]))
    }
    if(signal_num < ncol(data1)){data1[,(signal_num+1):feature_num] = rnorm(size*(feature_num-signal_num), 0, sqrt(sigma))}
  }else if(dist == "logit"){
    for(i in 1:length(ind_list)){
      data1[ind_list[[i]], 1:signal_num] = sapply(X = rnorm(length(ind_list[[i]])*signal_num,
                                                            eff_size[i],
                                               sqrt(sigma*variance_ratio[i])
                                              ),
                                         FUN = function(x) exp(x)/(1+exp(x))
                                        )
    }
    data1[,(signal_num+1):feature_num] = sapply(X = rnorm(size*(feature_num-signal_num), 0, sqrt(sigma)),
                             FUN = function(x) exp(x)/(1+exp(x)))

  }
  return(data1)
}




#n1 is # of correct weights for informative features
#n2 is # of correct weights for non-informative features
#r1, r2, uninfo_r1, uninfo_r2 define the magnitude of correct wights
get_weight = function(n1_ratio = 0.6,n2_ratio = 0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1,feature_num,signal_ratio){
  signal_num = feature_num*signal_ratio
  n1 = signal_num*n1_ratio
  n2 = (feature_num-signal_num)*n2_ratio
  ge_weight=abs(runif(feature_num,0,1))
  s=sample(1:signal_num,n1)
  ge_weight[s]=runif(n1,r1,r2)
  ss=sample((signal_num+1):feature_num,n2)
  ge_weight[ss]=runif(n2,uninfo_r1,uninfo_r2)
  return(ge_weight)
}

## 1.4 old functions ----
get_data_ori = function(size, index, data_divide, dist, eff_size,sigma,signal_ratio,feature_num,variance_ratio){
  data1 = matrix(0,size, feature_num)
  signal_num = feature_num*signal_ratio
  divide = lapply(X = str_split(str_split(data_divide,pattern = "/")[[1]], pattern = ""),
                  FUN = as.numeric
  )
  ind_list = list()
  for(i in 1: length(divide)){
    ind = NULL
    for(j in 1:length(divide[[i]])){
      ind = c(ind,index[[divide[[i]][[j]]]])
    }
    ind_list[[i]] = ind
  }
  if(length(eff_size) == 1){eff_size = rep(eff_size,length(divide))}
  if (dist == "normal"){
    for(i in 1:length(ind_list)){
      data1[ind_list[[i]], 1:signal_num] = rnorm(length(ind_list[[i]])*signal_num,
                                                 ((-1)^i)*ifelse(i<=2, eff_size[1], eff_size[1]*2),
                                                 sqrt(sigma*variance_ratio[i]))
    }
    data1[,(signal_num+1):feature_num] = rnorm(size*(feature_num-signal_num), 0, sqrt(sigma))
  }else if(dist == "logit"){
    for(i in 1:length(ind_list)){
      data1[ind_list[[i]], 1:signal_num] = sapply(X = rnorm(length(ind_list[[i]])*signal_num,
                                                            ((-1)^i)*ifelse(i<=2, eff_size[1], 2*eff_size[1]),
                                                            sqrt(sigma*variance_ratio[i])
      ),
      FUN = function(x) exp(x)/(1+exp(x))
      )
    }
    data1[,(signal_num+1):feature_num] = sapply(X = rnorm(size*(feature_num-signal_num), 0, sqrt(sigma)),
                                                FUN = function(x) exp(x)/(1+exp(x)))

  }
  return(data1)
}

















# 2 microbiome simulation ----
# 2.0 basic environment and packages
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
# 2.1 packages and main generation functions ----
# FUNCTION: microb_simulation
## input:
  # n: # samples
  # n_cluster: # clusters to be used
  # type: what type of info otu is used,
      # 1 means only 1 type of OTU (phylo_related), only works for different abundance level scenario(eta = NA, increased !=NA);
      # 3 means 3 types of OTU to be informative: 'Phy_related_abun', 'Phy_related_pres', 'Phy_unrelated_abun';
        # when generated using diff abundance level, no difference between abun and pres/abs
  # eta_threshold # threshold for affected
  # increased_perc
## output: data list with following elements
  # data
  # truelabel
  # dist_list




## generation main functions ----
microb_simulation = function(n = 100, # sample sizes for each cluster
                             n_cluster = 3, # number of clusters
                             type = 3, # type of info OTU
                             eta_perc = NA,
                             increased_perc = NA,
                             heter_B = F,
                             B = NA
){
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
    names(info_OTUs_ls_Phy_related_abun) = c('Phy_related_abun_1', 'Phy_related_abun_2', 'Phy_related_abun_3')

    # 1st set of info OTUs: the 3rd largest cluster. We will use abundance of these OTUs
    info_OTUs_ls_Phy_related_abun[[1]] = names(which(clustering == names(p_clus)[2]))
    info_OTUs_ls_Phy_related_abun[[2]] = names(which(clustering == names(p_clus)[3]))
    info_OTUs_ls_Phy_related_abun[[3]] = names(which(clustering == names(p_clus)[4]))

    info_OTUs_ls = info_OTUs_ls_Phy_related_abun

  }else if(type == 2){
    info_OTUs_ls = vector('list', length = 3) # list of info OTUs
    names(info_OTUs_ls) = c('Phy_related_pre/abs_1', 'Phy_related_pre/abs_2', 'Phy_related_pre/abs_3')

    # 1st set of info OTUs: the 3rd largest cluster. We will use abundance of these OTUs
    info_OTUs_ls[[1]] = names(which(clustering == names(p_clus)[18]))
    info_OTUs_ls[[2]] = names(which(clustering == names(p_clus)[19]))
    info_OTUs_ls[[3]] = names(which(clustering == names(p_clus)[20]))
    # length(info_OTUs_ls[[1]]) # 85 OTUs
    # sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1026
  }else if(type == 4){
    info_OTUs_ls = vector('list', length = 3)
    names(info_OTUs_ls) = c('Phy_related_abun', 'Phy_related_pres', 'Phy_unrelated_abun')

    info_OTUs_ls[[1]] = names(which(clustering == names(p_clus)[2]))
    info_OTUs_ls[[2]] = names(which(clustering == names(p_clus)[18]))
    temp_OTUs = names(sort(p_est, decreasing = T)[15:50])
    temp_TF = rep(T, length(temp_OTUs))
    for (j in 2:length(temp_OTUs)) { # to make sure OTUs are in different clusters
      if (clustering[temp_OTUs[j]] %in% clustering[temp_OTUs[1:(j-1)]]) { temp_TF[j] = F }
    }
    temp_OTUs = temp_OTUs[temp_TF]
    info_OTUs_ls[[3]] = temp_OTUs[-which(temp_OTUs %in% unlist(info_OTUs_ls[1:2]))]
  }

  # data generation ----
  if(!is.na(increased_perc)) {
    data = generate_diff_dirich_data(n = n,
                                     n_cluster = n_cluster,
                                     info_OTUs_ls = info_OTUs_ls,
                                     increased_perc = increased_perc)
  }else if (!is.na(eta_perc)){
    data = generate_case_control_data(n = n,
                                      n_cluster = n_cluster,
                                      info_OTUs_ls = info_OTUs_ls,
                                      eta_perc = eta_perc,
                                      type = type,
                                      heter_B = heter_B,
                                      B = B)
  }
  return(data)
}


generate_case_control_data = function(n = 100, # sample sizes for each cluster
                                      n_cluster = 3, # number of clusters
                                      # beta_group = 3, # effect size
                                      OTUnames = throat.tree$tip.label,
                                      info_OTUs_ls,
                                      eta_perc = 0.01, # percent of tail
                                      type = 3,
                                      heter_B = F, B = NA){
  for (k in 1:n_cluster) { # for each cluster
    info_OTUs = info_OTUs_ls[[k]] # informative OTUs
    ##### Generate data
    cum_OTUtab = NULL
    m = 1
    while (m < Inf) {
      OTUtab = sim_data_g_est(nSam = n, mu = 1000, size = 25, g_est = g_est)

      if(type == 1){# phylo-realted abun
        input_OTUtab = OTUtab / rowSums(OTUtab)
      }else if (type == 2){ # phylo-related pres/abs
        input_OTUtab = OTUtab
        input_OTUtab[input_OTUtab > 0] = 1
      }else if (type == 3){ # mixed type
        if (k == 2) {
            input_OTUtab = OTUtab
            input_OTUtab[input_OTUtab > 0] = 1 # change abundance to presence information
        } else {
          input_OTUtab = OTUtab / rowSums(OTUtab) # convert OTU counts to abundance
        }
      }

      if (k == 2) {
        if(type == 3){
          input_OTUtab = OTUtab
          input_OTUtab[input_OTUtab > 0] = 1 # change abundance to presence information
        }
      } else {
        input_OTUtab = OTUtab / rowSums(OTUtab) # convert OTU counts to abundance
      }
      X = input_OTUtab[, info_OTUs]
      if(heter_B){
        if(is.na(B)){B = runif(ncol(X),0,2)}else{
          B = B[colnames(input_OTUtab) %in% info_OTUs]
        }
      }else{
        B = rep(1, ncol(X))
      }
      eta = scale(rowSums(X%*%B))
      eta_threshold = quantile(eta,1-eta_perc, na.rm = T)
      temp_ids = which(eta >= eta_threshold) # samples with large eta
      Y = rep(0, length(eta)) # binary outcome
      Y[temp_ids] = rbinom(length(temp_ids), 1, 0.9) # if eta is large, then case with prob = 90%
      case_id = which(Y==1)
      if (is.null(cum_OTUtab) & length(case_id)>0) {
        cum_OTUtab = OTUtab[case_id, , drop = F] # cumulative OTU table
      } else {
        cum_OTUtab = rbind(cum_OTUtab, OTUtab[case_id, , drop = F])
      }
      if (nrow(cum_OTUtab) >= n) { break } else { m = m + 1 }
    }
    cum_OTUtab = cum_OTUtab[1:n, ]
    rownames(cum_OTUtab) = 1:n
    if (k == 1) {
      dataset = cbind(cum_OTUtab, rep(k, n))
      colnames(dataset) = c(OTUnames, 'Subtype')
    } else {
      dataset = rbind(dataset, cbind(cum_OTUtab, rep(k, n)))
    }
  }
  dataset = data.frame(dataset) # the dataset with OTU counts
  rownames(dataset) = 1:nrow(dataset)
  dist_list = dist_ls(dataset, OTUnames)
  list(data = dataset[,1:(ncol(dataset)-1)],
       truelabel = dataset$Subtype,
       dist_list = dist_list)
}


generate_diff_dirich_data=function(n = 100, # sample sizes for each cluster
                                   n_cluster = 3, # number of clusters
                                   beta_group = 3, # effect size,
                                   increased_perc = 0.1, # increased percentage on info OTUs
                                   OTUnames = throat.tree$tip.label,
                                   info_OTUs_ls){
  p_est_new = dd$pi # mle of prob
  g_est_new_ls = p_est_new_ls = list(p_est_new, p_est_new, p_est_new)
  for (ll in 1:3) { # for each info OTUs
    # IDs: 1: abundance increased set; 2~3: abundance decreased set
    if (ll == 1) { IDs = c(1,2,3) } else if (ll == 2) { IDs = c(2,3,1) } else { IDs = c(3,1,2)}
    current_p_est = p_est_new_ls[[ll]] # old proportion
    p_est_set1 = current_p_est[info_OTUs_ls[[IDs[1]]]] # old proportion of increased set
    p_est_set1 = p_est_set1 * (sum(p_est_set1) + increased_perc) / sum(p_est_set1) # adjust proportions
    p_est_set_23 = current_p_est[unlist(info_OTUs_ls[IDs[c(2,3)]])] # old proportion of decreased set
    p_est_set_23 = p_est_set_23 * ((sum(p_est_set_23)) - increased_perc) / sum(p_est_set_23) # adjust proportions
    ### Change values
    current_p_est[info_OTUs_ls[[IDs[1]]]] = p_est_set1 # record new proportions
    current_p_est[unlist(info_OTUs_ls[IDs[c(2,3)]])] = p_est_set_23 # record new proportions
    p_est_new_ls[[ll]] = current_p_est # record new proportions in a list
    g_est_new_ls[[ll]] = current_p_est * gplus # record new parameters in a list
    # sum(current_p_est) # new proportions; sum should be 1
  }
  for (k in 1:n_cluster) { # for each cluster
    OTUtab = sim_data_g_est(nSam = n, mu = 1000, size = 25, g_est = g_est_new_ls[[k]])
    if (k == 1) {
      dataset = cbind(OTUtab, rep(k, n))
      colnames(dataset) = c(OTUnames, 'Subtype')
    } else {
      dataset = rbind(dataset, cbind(OTUtab, rep(k, n)))
    }
  }
  dataset = data.frame(dataset) # the dataset with OTU counts
  rownames(dataset) = 1:nrow(dataset)
  dist_list = dist_ls(dataset, OTUnames)
  list(data = dataset[,1:(ncol(dataset)-1)],
       truelabel = dataset$Subtype,
       dist_list = dist_list)
}

generate_diff_dirich_data_single_set= function(n = 100, # sample sizes for each cluster
                                               n_cluster = 3, # number of clusters
                                               # beta_group = 3, # effect size,
                                               increased_perc = c(0.3,0.1,0), # increased percentage on info OTUs
                                               OTUnames = throat.tree$tip.label,
                                               info_OTUs_ls){
  p_est_new = dd$pi # mle of prob
  g_est_new_ls = p_est_new_ls = list(p_est_new, p_est_new, p_est_new)
  functional_set_ind = unlist(info_OTUs_ls)
  increasing_ind = functional_set_ind[1:floor(length(functional_set_ind)/3)]
  decreasing_ind = setdiff(functional_set_ind,increasing_ind)


  for(j in 1:length(increased_perc)){
    current_p_est = p_est_new_ls[[j]] # old proportion
    p_est_set_increase = current_p_est[increasing_ind] # old proportion of increased set
    p_est_set_increase = p_est_set_increase * (sum(p_est_set_increase) + increased_perc[j]) / sum(p_est_set_increase) # adjust proportions
    p_est_set_decrease = current_p_est[decreasing_ind] # old proportion of decreased set
    p_est_set_decrease = p_est_set_decrease * ((sum(p_est_set_decrease)) - increased_perc[j]) / sum(p_est_set_decrease) # adjust proportions
    ### Change values
    current_p_est[increasing_ind] = p_est_set_increase # record new proportions
    current_p_est[decreasing_ind] = p_est_set_decrease # record new proportions
    p_est_new_ls[[j]] = current_p_est # record new proportions in a list
    g_est_new_ls[[j]] = current_p_est * gplus # record new parameters in a list
    cat(j,"sum of p is ", sum(current_p_est),"\n") # new proportions; sum should be 1)
  }

  for (k in 1:n_cluster) { # for each cluster
    OTUtab = sim_data_g_est(nSam = n, mu = 1000, size = 25, g_est = g_est_new_ls[[k]])
    if (k == 1) {
      dataset = cbind(OTUtab, rep(k, n))
      colnames(dataset) = c(OTUnames, 'Subtype')
    } else {
      dataset = rbind(dataset, cbind(OTUtab, rep(k, n)))
    }
  }
  dataset = data.frame(dataset) # the dataset with OTU counts
  rownames(dataset) = 1:nrow(dataset)
  dist_list = dist_ls(dataset, OTUnames)
  list(data = dataset[,1:(ncol(dataset)-1)],
       truelabel = dataset$Subtype,
       dist_list = dist_list)
}






# 3 2 moons simulation ----
simu_2moon = function (r = 1, n = 200, sigma = 0.1, levels = NULL,
                       seed = NULL)
{
  if(!is.null(seed)){
    set.seed(seed)
  }

  if (length(n) == 1)
    n = rep(n, 2)
  if (is.null(levels))
    levels = paste("Class", 1:2)
  alpha1 = stats::runif(n[1], 0, pi)
  alpha2 = stats::runif(n[2], 0, pi)
  k = c(rep(1, n[1]), rep(2, n[2]))
  sx = r/2
  sy = r/6
  x = c(r * cos(alpha1) + sx, r * cos(alpha2) - sx)
  y = c(r * sin(alpha1) - sy, r * -sin(alpha2) + sy)
  d = cbind.data.frame(x, y, factor(k, labels = levels))
  noise = matrix(stats::rnorm(2 * sum(n), 0, sigma), ncol = 2)
  d[, 1:2] = d[, 1:2] + noise
  colnames(d) = c("X", "Y", "Class")
  return(d)
}

## utils ----
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


