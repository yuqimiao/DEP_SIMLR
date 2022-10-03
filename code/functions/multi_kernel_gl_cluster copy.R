
multi_kernels_gl = function (x = NA,
                             Diff = NA,
                             sigma =seq(2, 1, -0.25),
                             allk = seq(10,30,2),
                             cores.ratio = 0.5,
                             is.square = F,
                             standardize = 0,
                             feature_weight = NA) {


  N = dim(x)[1]
  KK = 0
  mu = seq(0.1,1, 0.1)
  allk = allk

  # distance calc and sort
  if(is.na(Diff)){
    if(!is.na(feature_weight)){
      if(is.square){
        Diff = dist2_w(x,x, feature_weight)
      }else{
        Diff = dist2_w(x,x, feature_weight)^(1/2)
      }
    }else{
      if(is.square){
        Diff = dist2(x,x)
      }else{
        Diff = dist2(x,x)^(1/2)
      }
    }
  }



  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    .libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
    library(Matrix)
    library(abSNF)
    library(SNFtool)
    D_kernels = function(Diff_fun = NA,
                         Diff_sort_fun = NA,
                         allk_fun = 30,
                         sigma_fun = 2) {

      # calc symmetric mu_{ij} for every subject
      TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      ## calc mean distance for the first k neighbors of every subjects(row), length = k
      TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
      Sig = Sig + t(Sig)
      Sig = Sig / 2
      Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      ## make the symmetric mu_{ij}

      W = dnorm(Diff_fun,0,sigma_fun*Sig)
      # Output the symmetric kernel by matrix
      W_ori = D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)

      # if is.SP, using dij = kii+kjj-2kijto calculate distance based on kernel value
      # SIMLR proceed
      K = D_Kernels
      k = 1/sqrt(diag(K)+1)
      # diagnal of the matrix is the highest similarity
      G = K * (k %*% t(k))
      ## with diag(K) all ~ 0, k is just n 1s, thus G is K
      G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
      G2 = t(G1)
      # Use the average difference btw 2 self similarities and pairwise similarity as kenels
      D_Kernels_tmp = (G1 + G2 - 2*G)/2
      D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
      D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
      Kernels = max(as.matrix(D_Kernels))-D_Kernels
      return(list(ori_Kernels = W_ori, Kernels = Kernels, D_Kernels = D_Kernels))
    }
  })

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l,Diff_fun = Diff,
                                                                  Diff_sort_fun = Diff_sort,
                                                                  B = allk,
                                                                  sig_df_list = sigma,
                                                                  mu_aff_list = mu,
                                                                  KK_fun = KK) {
    for (j in 1:length(sig_df_list)){
      name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
      D_Kernels[[name]] = Matrix(D_kernels(Diff_fun = Diff_fun,
                                           Diff_sort_fun = Diff_sort_fun,
                                           allk_fun = B[l],
                                           sigma_fun = sig_df_list[[j]])[[1]])
    }

    # for(t in 1:length(mu_aff_list)){
    #   name = paste("aff_B",B[l], "_mu",mu_aff_list[[t]],sep = "" )
    #   D_Kernels[[name]] = Matrix(affinityMatrix(dist2(x,x)^(1/2),K = B[l], sigma = mu_aff_list[[t]]))
    # }


    return(D_Kernels)
  }
  ))
  stopCluster(cl)
  return(D_Kernels)
}


# old version ----
# multi_kernels_gl = function (x, sigma =seq(2, 1, -0.25),
#                              allk = seq(10,30,2),
#                              cores.ratio = 0.5,
#                              is.square = F,
#                              standardize = 0,
#                              feature_weight = NA) {
#
#
#   N = dim(x)[1]
#   KK = 0
#   mu = seq(0.1,1, 0.1)
#   allk = allk
#   cores = as.integer(cores.ratio * (detectCores() - 1))
#   if (cores < 1 || is.na(cores) || is.null(cores)) {
#     cores = 1
#   }
#   cl = makeCluster(cores)
#   clusterEvalQ(cl, {
#     .libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
#     library(Matrix)
#     library(abSNF)
#     library(SNFtool)
#     D_kernels = function(x_fun = NA,
#                          Diff_fun = NA,
#                          allk_fun = 30,
#                          sigma_fun = 2,
#                          is.square = F,
#                          is.normal = T,
#                          weight = NA,
#                          is.SP = T,
#                          standardization = standardize) {
#       if(is.na(Diff_fun)){
#         if(standardization == 1){
#           x_fun = standardNormalization(x_fun)
#         }else if(standardization == 2){
#           x_fun = apply(x_fun, 2, function(y) {(y-min(y))/(max(y)-min(y))})
#         }
#
#
#         if(allk_fun<(nrow(x_fun))) {
#           # distance calc and sort
#           if(!is.na(weight)){
#             if(is.square){
#               Diff_fun = dist2_w(x_fun,x_fun, weight)^2
#             }else{
#               Diff_fun = dist2_w(x_fun,x_fun, weight)^(1/2)
#             }
#           }else{
#             if(is.square){
#               Diff_fun = dist2(x_fun,x_fun)^2
#             }else{
#               Diff_fun = dist2(x_fun,x_fun)^(1/2)
#             }
#           }
#         }
#       }
#       diag(Diff_fun) = 0
#       Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort))
#       ## sort every column and transpose => ith row is the sorted distance btw ith subject and others
#
#       # calc symmetric mu_{ij} for every subject
#       TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
#       ## calc mean distance for the first k neighbors of every subjects(row), length = k
#       TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
#       Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
#       Sig = Sig + t(Sig)
#       Sig = Sig / 2
#       Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
#       Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
#       Sig = Sig * Sig_valid + .Machine$double.eps
#       ## make the symmetric mu_{ij}
#
#       if(is.normal){
#         W = dnorm(Diff_fun,0,sigma_fun*Sig)
#       }else{
#         W = dnorm(Diff_fun,0,sigma_fun*Sig)*sigma_fun*Sig
#       }
#
#       # Output the symmetric kernel by matrix
#       W_ori = D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
#
#       # if is.SP, using dij = kii+kjj-2kijto calculate distance based on kernel value
#       if(is.SP){
#         # SIMLR proceed
#         #update from SIMLR function (??)
#         K = D_Kernels
#         k = 1/sqrt(diag(K)+1)
#         # diagnal of the matrix is the highest similarity
#         G = K * (k %*% t(k))
#         ## with diag(K) all ~ 0, k is just n 1s, thus G is K
#         G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
#         G2 = t(G1)
#         # Use the average difference btw 2 self similarities and pairwise similarity as kenels
#         D_Kernels_tmp = (G1 + G2 - 2*G)/2
#         D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
#         D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
#       }
#       Kernels = max(as.matrix(D_Kernels))-D_Kernels
#       return(list(ori_Kernels = W_ori, Kernels = Kernels, D_Kernels = D_Kernels))
#     }
#   })
#
#   D_Kernels = list()
#   D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, x_fun = x, B = allk,
#                                                                   sig_df_list = sigma, mu_aff_list = mu, KK_fun = KK,
#                                                                   is.square_fun = is.square, weight_fun = feature_weight) {
#     if (B[l] < (nrow(x_fun) - 1)) {
#       for (j in 1:length(sig_df_list)){
#         name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
#         D_Kernels[[name]] = Matrix(D_kernels(x_fun = x,
#                                              allk_fun = B[l],
#                                              sigma_fun = sig_df_list[[j]],
#                                              standardization = 0,
#                                              is.square = is.square_fun,
#                                              weight = weight_fun)[[1]])
#       }
#
#       # for(t in 1:length(mu_aff_list)){
#       #   name = paste("aff_B",B[l], "_mu",mu_aff_list[[t]],sep = "" )
#       #   D_Kernels[[name]] = Matrix(affinityMatrix(dist2(x,x)^(1/2),K = B[l], sigma = mu_aff_list[[t]]))
#       # }
#
#     }
#
#     return(D_Kernels)
#   }
#   ))
#   stopCluster(cl)
#   return(D_Kernels)
# }

## for large feature size ----
multi_kernels_gl_t = function (x = NA, dist = NA, sigma =seq(2, 1, -0.1), cores.ratio = 0.5,is.square = F, standardize = 0) {
  KK = 0
  mu = seq(0.1,1,0.1)
  allk = c(10,30)

  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }

  if(is.na(dist)){
    N = dim(x)[2]
    Diff = dist2_t(x)^(1/2)
  }else{
    N = dim(dist)[1]
    Diff = dist^(1/2)
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    .libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
    library(Matrix)
    library(SNFtool)
    D_kernels = function(x_fun = NA,
                         Diff_fun = NA,
                         allk_fun = 30,
                         sigma_fun = 2,
                         is.square = F,
                         is.normal = T,
                         is.SP = T,
                         standardization = standardize) {
      if(is.na(Diff_fun)){
        if(standardization == 1){
          x_fun = standardNormalization(x_fun)
        }else if(standardization == 2){
          x_fun = apply(x_fun, 2, function(y) {(y-min(y))/(max(y)-min(y))})
        }


        if(allk_fun<(nrow(x_fun))) {
          # distance calc and sort
          if(is.square){
            Diff_fun = dist2_t(x_fun)^2
          }else{
            Diff_fun = dist2_t(x_fun)^(1/2)
          }
        }
      }
      diag(Diff_fun) = 0
      Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort))
      ## sort every column and transpose => ith row is the sorted distance btw ith subject and others

      # calc symmetric mu_{ij} for every subject
      TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      ## calc mean distance for the first k neighbors of every subjects(row), length = k
      TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
      Sig = Sig + t(Sig)
      Sig = Sig / 2
      Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      ## make the symmetric mu_{ij}

      if(is.normal){
        W = dnorm(Diff_fun,0,sigma_fun*Sig)
      }else{
        W = dnorm(Diff_fun,0,sigma_fun*Sig)*sigma_fun*Sig
      }

      # Output the symmetric kernel by matrix
      W_ori = D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)

      # if is.SP, using dij = kii+kjj-2kijto calculate distance based on kernel value
      if(is.SP){
        # SIMLR proceed
        #update from SIMLR function (??)
        K = D_Kernels
        k = 1/sqrt(diag(K)+1)
        # diagnal of the matrix is the highest similarity
        G = K * (k %*% t(k))
        ## with diag(K) all ~ 0, k is just n 1s, thus G is K
        G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
        G2 = t(G1)
        # Use the average difference btw 2 self similarities and pairwise similarity as kenels
        D_Kernels_tmp = (G1 + G2 - 2*G)/2
        D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
        D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
      }
      Kernels = max(as.matrix(D_Kernels))-D_Kernels
      return(list(ori_Kernels = W_ori, Kernels = Kernels, D_Kernels = D_Kernels))
    }

  })

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, Diff_fun = Diff, B = allk,
                                                                  sig_df_list = sigma, mu_aff_list = mu, KK_fun = KK,
                                                                  is.square_fun = is.square) {
    if (B[l] < (nrow(Diff_fun) - 1)) {
      for (j in 1:length(sig_df_list)){
        name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
        D_Kernels[[name]] = Matrix(D_kernels(Diff_fun = Diff_fun,
                                             allk_fun = B[l],
                                             sigma_fun = sig_df_list[[j]],
                                             standardization = 0,
                                             is.square = is.square_fun)[[1]])
      }

      # for(t in 1:length(mu_aff_list)){
      #   name = paste("aff_B",B[l], "_mu",mu_aff_list[[t]],sep = "" )
      #   D_Kernels[[name]] = Matrix(affinityMatrix(dist2(x,x)^(1/2),K = B[l], sigma = mu_aff_list[[t]]))
      # }

    }

    return(D_Kernels)
  }
  ))
  stopCluster(cl)
  return(D_Kernels)
}







