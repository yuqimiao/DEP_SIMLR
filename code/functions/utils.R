# Core utils ----

"dist2" = function( x, c = NA ) {

  # set the parameters for x
  if(is.na(c)) {
    c = x
  }

  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }

  # compute the distance
  dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) +
    (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) -
    2 * (x%*%t(c))

  dist[dist<0] = 0

  return(dist)

}
standardNormalization = function (x){
  x = as.matrix(x)
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd == 0] = 1
  xNorm = t((t(x) - mean)/sd)
  return(xNorm)
}
multi_kernels_gl = function (x = NA,
                             Diff = NA,
                             sigma =seq(2, 1, -0.25),
                             allk = seq(10,30,2),
                             cores.ratio = 0.5,
                             is.square = F,
                             with_coef = T,
                             standardize = 0,
                             feature_weight = NA) {


  # N = dim(x)[1]
  KK = 0
  mu = seq(0.1,1, 0.1)
  allk = allk

  # distance calc and sort
  if(is.na(Diff)){
    if(!is.na(feature_weight)){
      Diff = dist2_w(x,x, feature_weight)

    }else{
      Diff = dist2(x,x)
    }
  }

  Diff[Diff<0] = 0
  if(!is.square){
    Diff = Diff^(1/2)
  }

  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    library(Matrix)
    library(abSNF)
    library(SNFtool)
    D_kernels = function(Diff_fun = NA,
                         Diff_sort_fun = NA,
                         allk_fun = 30,
                         sigma_fun = 2,
                         with_coef_fun = T) {

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

      if(with_coef_fun){
        W = dnorm(Diff_fun,0,sigma_fun*Sig)
      }else{
        W = dnorm(Diff_fun,0,sigma_fun*Sig)*Sig
      }
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
                                                                  KK_fun = KK,
                                                                  with_coef_fun = with_coef) {
    for (j in 1:length(sig_df_list)){
      name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
      D_Kernels[[name]] = Matrix(D_kernels(Diff_fun = Diff_fun,
                                           Diff_sort_fun = Diff_sort_fun,
                                           allk_fun = B[l],
                                           sigma_fun = sig_df_list[[j]],
                                           with_coef_fun = with_coef_fun)[[1]])
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

multi_kernels_2 = function (x = NA,
                           Diff = NA,
                           sigma = sqrt(seq(2, 1, -0.25)),
                           allk = seq(10,30,2),
                           cores.ratio = 0.5,
                           is.square = F,
                           # with_coef = T,
                           standardize = 0,
                           feature_weight = NA,
                           kernel_f = "SIMLR") {


  # N = dim(x)[1]
  KK = 0
  mu = seq(0.1,1, 0.1)
  allk = allk

  # distance calc and sort
  if(is.na(Diff)){
    if(!is.na(feature_weight)){
      Diff = dist2_w(x,x, feature_weight)

    }else{
      Diff = dist2(x,x)
    }
  }

  Diff[Diff<0] = 0
  if(!is.square){
    Diff = Diff^(1/2)
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
    D_kernels_2 = function(Diff_fun = NA,
                           Diff_sort_fun = NA,
                           allk_fun = 30,
                           sigma_fun = sqrt(2),
                           # with_coef_fun = T,
                           kernel_form = "SIMLR"# options: SNF, SNF_2, SIMLR, SIMLR_2, uni_denorm
    ) {

      # calc symmetric mu_{ij} for every subject
      TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      ## calc mean distance for the first k neighbors of every subjects(row), length = k
      TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
      Sig[Sig<0] = 0
      if(kernel_form %in% c("SIMLR", "SIMLR_2")){
        Sig = Sig + t(Sig)
        Sig = Sig / 2
      }else if(kernel_form %in% c("SNF", "SNF_2")){
        Sig = (Sig+t(Sig)+Diff_fun)/3
        if(kernel_form == "SNF"){
          Sig = sqrt(Sig)
        }
      }else if(kernel_form == "uni_denorm"){
        ns = nrow(Diff_fun)
        Sig = matrix(mean(Diff_fun), ns, ns)
      }
      Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      ## make the symmetric mu_{ij}

      if(kernel_form %in% c("SNF","SNF_2")){
        sigma_fun = sqrt(sigma_fun/2)
      }

      if(kernel_form %in% c("SIMLR","SNF_2" )){
        W = dnorm(Diff_fun,0,sigma_fun*Sig)
      }else{
        W = dnorm(Diff_fun,0,sigma_fun*Sig)*sigma_fun*Sig*sqrt(2*pi)
      }

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
                                                                  kernel_form = kernel_f) {
    for (j in 1:length(sig_df_list)){
      name = paste("dk_B",B[l], "_sig",sig_df_list[j],sep = "" )
      D_Kernels[[name]] = Matrix(D_kernels_2(Diff_fun = Diff_fun,
                                             Diff_sort_fun = Diff_sort_fun,
                                             allk_fun = B[l],
                                             sigma_fun = sig_df_list[j],
                                             kernel_form = kernel_form)[[1]])
    }



    return(D_Kernels)
  }
  ))
  stopCluster(cl)
  return(D_Kernels)
}
# old version with more flexible square,standard setting
multi_kernels_gl_more_par = function (x, sigma =seq(2, 1, -0.1),allk = c(10,30), cores.ratio = 0.5,is.square = F, standardize = 0, feature_weight = NA) {


  N = dim(x)[1]
  KK = 0
  mu = seq(0.1,1, 0.1)
  allk = allk
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    library(Matrix)
    library(abSNF)
    library(SNFtool)
    D_kernels = function(x_fun = NA,
                         Diff_fun = NA,
                         allk_fun = 30,
                         sigma_fun = 2,
                         is.square = F,
                         is.normal = T,
                         weight = NA,
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
          if(!is.na(weight)){
            if(is.square){
              Diff_fun = dist2_w(x_fun,x_fun, weight)
            }else{
              Diff_fun = dist2_w(x_fun,x_fun, weight)^(1/2)
            }
          }else{
            if(is.square){
              Diff_fun = dist2(x_fun,x_fun)
            }else{
              Diff_fun = dist2(x_fun,x_fun)^(1/2)
            }
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
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, x_fun = x, B = allk,
                                                                  sig_df_list = sigma, mu_aff_list = mu, KK_fun = KK,
                                                                  is.square_fun = is.square, weight_fun = feature_weight) {
    if (B[l] < (nrow(x_fun) - 1)) {
      for (j in 1:length(sig_df_list)){
        name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
        D_Kernels[[name]] = Matrix(D_kernels(x_fun = x,
                                             allk_fun = B[l],
                                             sigma_fun = sig_df_list[[j]],
                                             standardization = 0,
                                             is.square = is.square_fun,
                                             weight = weight_fun)[[1]])
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

## for large feature size
"dist2_t" = function(x) { ## x should have features at row and samples at column

  # compute the dimension
  n1 = ncol(x)
  d1 = nrow(x)
  x = as.matrix(x)
  # compute the distance
  square_sum_mat = matrix(1,nrow=n1) %*% matrix(apply(x^2,MARGIN=2,FUN=sum),nrow = 1)
  dist = square_sum_mat+t(square_sum_mat) - 2 * (t(x)%*%x)
  diag(dist) = 0
  return(dist)

}
standardNormalization_real = function(x){
  direction = ifelse(dim(x)[1]>=dim(x)[2],1,2)
  mean = apply(x,1,mean)
  sd = apply(x,1,sd)
  (x-mean)/sd
}
multi_kernels_gl_t = function (x = NA, dist = NA, sigma =seq(2, 1, -0.1), cores.ratio = 0.5,is.square = F, standardize = 1) {
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




# using the std of data: (y-min(y))/(max(y)-min(y)) make sure entries in (0,1)
# using the defination as formula in SIMLR
D_kernels = function(x_fun = NA,
                     Diff_fun = NA,
                     allk_fun = 30,
                     sigma_fun = 2,
                     is.square = F,
                     is.normal = T,
                     is.SP = T,
                     standardization = 1) {
  if(is.na(Diff_fun)){
    if(standardization == 1){
      x_fun = standardNormalization(x_fun)
    }else if(standardization == 2){
      x_fun = apply(x_fun, 2, function(y) {(y-min(y))/(max(y)-min(y))})
    }


    if(allk_fun<(nrow(x_fun))) {
      # distance calc and sort
      if(is.square){
        Diff_fun = dist2(x_fun,x_fun)^2
      }else{
        Diff_fun = dist2(x_fun,x_fun)^(1/2)
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

dominateset <- function(xx,KK=20) {
  ### This function outputs the top KK neighbors.
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
  }

  return(A)
}

kernel_similarity = function(kernel_list){
  term = rep(0, length(kernel_list))
  for(i in 1:length(kernel_list)){
    for(j in 1:length(kernel_list)){
      term[i] = term[i]+sum(kernel_list[[i]]*kernel_list[[j]])
    }
  }
  return(term)
}

"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {
  ## Calculate eigen vector and eigen value given dimension
  # set the needed parameters
  if(is.na(c)) {
    c = dim(A)[1]
  }
  if(c>dim(A)[1]) {
    c = dim(A)[1]
  }
  if(is.na(isMax)) {
    isMax = 1
  }
  if(is.na(isSym)) {
    isSym = 1
  }

  # compute the eigenvalues and eigenvectors of A
  if(isSym==1) {
    eigen_A = eigen(A,symmetric=TRUE)
  }
  else {
    eigen_A = eigen(A)
  }
  v = eigen_A$vectors
  d = eigen_A$values

  # sort the eigenvectors
  if(isMax == 0) {
    eigen_A_sorted = sort(d,index.return=TRUE)
  }
  else {
    eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
  }
  d1 = eigen_A_sorted$x
  idx = eigen_A_sorted$ix
  idx1 = idx[1:c]

  # compute the results
  eigval = d[idx1]
  eigvec = Re(v[,idx1])
  eigval_full = d[idx]

  return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))

}

dist_kernels = function(kernel){
  K = kernel
  k = 1/sqrt(diag(K)+1)
  ## diagnal of the matrix is the highest similarity
  G = K * (k %*% t(k))
  ## with diag(K) all ~ 0, k is just n 1s, thus G is K
  G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
  G2 = t(G1)
  ## Use the average difference btw 2 self similarities and pairwise similarity as kenels
  D_Kernels_tmp = (G1 + G2 - 2*G)/2
  D_Kernels = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
  D_Kernels = (as.matrix(D_Kernels)+t(as.matrix(D_Kernels)))/2
  return(D_Kernels)
}

silhouette_similarity = function (group, similarity_matrix) {
  similarity_matrix = as.matrix(similarity_matrix)
  similarity_matrix <- (similarity_matrix + t(similarity_matrix))/2
  diag(similarity_matrix) = 0
  normalize <- function(X) X/rowSums(X)
  similarity_matrix <- normalize(similarity_matrix)
  n <- length(group)
  if (!all(group == round(group)))
    stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if (k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if (doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  wds <- matrix(NA, n, 3, dimnames = list(names(group), c("cluster",
                                                          "neighbor", "sil_width")))
  for (j in 1:k) {
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index,
                                            drop = FALSE], 2, function(r) tapply(r, group[!index],
                                                                                 mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index, "neighbor"] <- cluster_id[-j][maxC]
    s.i <- if (Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj -
                                                         1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i)/pmax(b.i, a.i), 0)
    }
    else 0
    wds[index, "sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}
# tunning functions ----
# Input
  # * scen_kernel_list: kernel list for certain sceanrio
  # * par_list: parameter list to be  used,
  #   * every element in this list should be a list with 5 components
  #     "rho","alpha0","alpha","beta","gamma"
# output:
  #* res: the result collection for each paramter combination
  #* par_eva: parameter evaluation using nmi, rand and silhouette.
dfl_simlr_tune = function(scen_kernel_list, par_list, truelabel, c = 4,
                          new_beta = F,
                          is.square = F,
                          cluster = "spectral",
                          network.diffusion = F,
                          normalize = 0,
                          k = 30){
  res = list()
  par_eva = NULL
  for(i in 1:length(par_list)){
    print(i)
    par_name = paste("alpha0.",par_list[[i]]$alpha0,
                     "_alpha.",par_list[[i]]$alpha,
                     "_rho.",par_list[[i]]$rho,
                     "_beta.",par_list[[i]]$beta,
                     "_gamma.",par_list[[i]]$gamma,
                     sep = "")
    tryCatch({
      res[[par_name]] = dfl_simlr(kernel_list = scen_kernel_list,
                                  c=c,
                                  k = k,
                                  alpha0 = par_list[[i]]$alpha0,
                                  alpha = par_list[[i]]$alpha,
                                  rho = par_list[[i]]$rho,
                                  beta = par_list[[i]]$beta,
                                  gamma = par_list[[i]]$gamma,
                                  normalize = normalize,
                                  new_beta = new_beta,
                                  network.diffusion = network.diffusion)
    },
    error = function(e) {
      cat(par_name," can't run, with error:", conditionMessage(e))
      })

    if(!is.null(res[[par_name]])){
      if(cluster == "spectral"){
        cluster_scen1_tmp = spectralClustering(as.matrix(res[[par_name]]$S), K = c)
      }else{
        cluster_scen1_tmp = kmeans(as.matrix(res[[par_name]]$L), centers = c, nstart = 200)$cluster
      }
      if(class(truelabel) == "list"){
        truelabel_i = truelabel[[i]]
      }else{
        truelabel_i = truelabel
      }

      nmi_tmp=igraph::compare(truelabel_i, cluster_scen1_tmp, method = "nmi")
      rand_tmp=igraph::compare(truelabel_i, cluster_scen1_tmp, method = "rand")
      silhouette_tmp = summary(silhouette_similarity(cluster_scen1_tmp, res[[par_name]]$S))$avg.width
      par_eva = rbind(par_eva, c(par_name = par_name,
                                 unlist(par_list[[i]]),
                                 nmi = nmi_tmp,
                                 rand = rand_tmp,
                                 silhouette = silhouette_tmp))
    }else{
      par_eva = rbind(par_eva, c(par_name = par_name,
                                 unlist(par_list[[i]]),
                                 nmi = "NULL",
                                 rand = "NULL",
                                 silhouette = "NULL"))
    }
  }

  return(list(res = res,par_eva = par_eva))
}


# CIMLR functions  ----
CIMLR_dist = function(dist_list,c, no.dim = NA, k = 10, cores.ratio = 1)
{
  if (is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  NITER = 30
  num = ncol(dist_list[[1]])
  r = -1
  beta = 0.8
  cat("Computing the multiple Kernels.\n")
  D_Kernels = unlist(lapply(dist_list, multiple.kernel.cimlr.dist))
  alphaK = 1/rep(length(D_Kernels), length(D_Kernels))
  distX = array(0, c(dim(D_Kernels[[1]])[1], dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX/length(D_Kernels)
  res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x,
                                                               index.return = TRUE)))
  distX1 = array(0, c(nrow(distX), ncol(distX)))
  idx = array(0, c(nrow(distX), ncol(distX)))
  for (i in 1:nrow(distX)) {
    distX1[i, ] = res[[i]]$x
    idx[i, ] = res[[i]]$ix
  }
  A = array(0, c(num, num))
  di = distX1[, 2:(k + 2)]
  rr = 0.5 * (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1,
                                      FUN = sum))
  id = idx[, 2:(k + 2)]
  numerator = (apply(array(0, c(length(di[, k + 1]), dim(di)[2])),
                     MARGIN = 2, FUN = function(x) {
                       x = di[, k + 1]
                     }) - di)
  temp = (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1, FUN = sum) +
            .Machine$double.eps)
  denominator = apply(array(0, c(length(temp), dim(di)[2])),
                      MARGIN = 2, FUN = function(x) {
                        x = temp
                      })
  temp = numerator/denominator
  a = apply(array(0, c(length(t(1:num)), dim(di)[2])), MARGIN = 2,
            FUN = function(x) {
              x = 1:num
            })
  A[cbind(as.vector(a), as.vector(id))] = as.vector(temp)
  if (r <= 0) {
    r = mean(rr)
  }
  lambda = max(mean(rr), 0)
  A[is.nan(A)] = 0
  S0 = max(max(distX)) - distX
  cat("Performing network diffusion.\n")
  S0 = network.diffusion(S0, k)
  S0 = dn.cimlr(S0, "ave")
  S = (S0 + t(S0))/2
  D0 = diag(apply(S, MARGIN = 2, FUN = sum))
  L0 = D0 - S
  eig1_res = eig1(L0, c, 0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  F_eig1 = dn.cimlr(F_eig1, "ave")
  converge = vector()
  for (iter in 1:NITER) {
    cat("Iteration: ", iter, "\n")
    distf = L2_distance_1(t(F_eig1), t(F_eig1))
    A = array(0, c(num, num))
    b = idx[, 2:dim(idx)[2]]
    a = apply(array(0, c(num, ncol(b))), MARGIN = 2, FUN = function(x) {
      x = 1:num
    })
    inda = cbind(as.vector(a), as.vector(b))
    ad = (distX[inda] + lambda * distf[inda])/2/r
    dim(ad) = c(num, ncol(b))
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx", c_input, c_output))
    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    S = (1 - beta) * A + beta * S
    S = network.diffusion(S, k)
    S = (S + t(S))/2
    D = diag(apply(S, MARGIN = 2, FUN = sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L, c, 0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    F_eig1 = dn.cimlr(F_eig1, "ave")
    F_eig1 = (1 - beta) * F_old + beta * F_eig1
    evs_eig1 = cbind(evs_eig1, ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps + D_Kernels[[i]]) * (S +
                                                         .Machine$double.eps)
      DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
    }
    alphaK0 = umkl.cimlr2(DD)
    alphaK0 = alphaK0/sum(alphaK0)
    alphaK = (1 - beta) * alphaK + beta * alphaK0
    alphaK = alphaK/sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c + 1)])
    converge[iter] = fn2 - fn1
    if (iter < 10) {
      if (ev_eig1[length(ev_eig1)] > 1e-06) {
        lambda = 1.5 * lambda
        r = r/1.01
      }
    }
    else {
      if (converge[iter] > 1.01 * converge[iter - 1]) {
        S = S_old
        if (converge[iter - 1] > 0.2) {
          warning("Maybe you should set a larger value of c.")
        }
        break
      }
    }
    S_old = S
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }
    res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x,
                                                                 index.return = TRUE)))
    distX1 = array(0, c(nrow(distX), ncol(distX)))
    idx = array(0, c(nrow(distX), ncol(distX)))
    for (i in 1:nrow(distX)) {
      distX1[i, ] = res[[i]]$x
      idx[i, ] = res[[i]]$ix
    }
  }
  LF = F_eig1
  D = diag(apply(S, MARGIN = 2, FUN = sum))
  L = D - S
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  if (length(no.dim) == 1) {
    U_index = seq(ncol(U), (ncol(U) - no.dim + 1))
    F_last = tsne(S, k = no.dim, initial_config = U[, U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U), (ncol(U) - no.dim[i] + 1))
      F_last[i] = list(tsne(S, k = no.dim[i], initial_config = U[,
                                                                 U_index]))
    }
  }
  execution.time = proc.time() - ptm
  cat("Performing Kmeans.\n")
  y = kmeans(F_last, c, nstart = 200)
  ydata = tsne(S)
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF
  return(results)
}

multiple.kernel.cimlr.dist = function(dist, cores.ratio = 1 ) {


  # compute some parameters from the kernels
  N = dim(dist)[1]
  KK = 0
  sigma = seq(2,1,-0.25)

  # compute and sort Diff
  Diff = dist
  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = seq(10,30,2)

  # setup a parallelized estimation of the kernels
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }

  cl = makeCluster(cores)

  clusterEvalQ(cl, {library(Matrix)})

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl,1:length(allk),fun=function(l,Diff_sort_fun=Diff_sort,allk_fun=allk,
                                                              Diff_fun=Diff,sigma_fun=sigma,KK_fun=KK) {
    if(allk_fun[l]<(nrow(Diff)-1)) {
      TT = apply(Diff_sort_fun[,2:(allk_fun[l]+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      TT = matrix(data = TT, nrow = length(TT), ncol = 1)
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
      Sig = Sig + t(Sig)
      Sig = Sig / 2
      Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      for (j in 1:length(sigma_fun)) {
        W = dnorm(Diff_fun,0,sigma_fun[j]*Sig)
        D_Kernels[[KK_fun+l+j]] = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
      }
      return(D_Kernels)
    }
  }))

  stopCluster(cl)

  # compute D_Kernels
  for (i in 1:length(D_Kernels)) {
    K = D_Kernels[[i]]
    k = 1/sqrt(diag(K)+1)
    G = K * (k %*% t(k))
    G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
    G2 = t(G1)
    D_Kernels_tmp = (G1 + G2 - 2*G)/2
    D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    D_Kernels[[i]] = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
  }

  return(D_Kernels)

}

# compute and returns the multiple kernel
"multiple.kernel.cimlr" = function( x, cores.ratio = 1 ) {

  # set the parameters
  kernel.type = list()
  kernel.type[1] = list("poly")
  kernel.params = list()
  kernel.params[1] = list(0)

  # compute some parameters from the kernels
  N = dim(x)[1]
  KK = 0
  sigma = seq(2,1,-0.25)

  # compute and sort Diff
  Diff = dist2.cimlr(x)
  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = seq(10,30,2)

  # setup a parallelized estimation of the kernels
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }

  cl = makeCluster(cores)

  clusterEvalQ(cl, {library(Matrix)})

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl,1:length(allk),fun=function(l,x_fun=x,Diff_sort_fun=Diff_sort,allk_fun=allk,
                                                              Diff_fun=Diff,sigma_fun=sigma,KK_fun=KK) {
    if(allk_fun[l]<(nrow(x_fun)-1)) {
      TT = apply(Diff_sort_fun[,2:(allk_fun[l]+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      TT = matrix(data = TT, nrow = length(TT), ncol = 1)
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
      Sig = Sig + t(Sig)
      Sig = Sig / 2
      Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      for (j in 1:length(sigma_fun)) {
        W = dnorm(Diff_fun,0,sigma_fun[j]*Sig)
        D_Kernels[[KK_fun+l+j]] = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
      }
      return(D_Kernels)
    }
  }))

  stopCluster(cl)

  # compute D_Kernels
  for (i in 1:length(D_Kernels)) {
    K = D_Kernels[[i]]
    k = 1/sqrt(diag(K)+1)
    G = K * (k %*% t(k))
    G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
    G2 = t(G1)
    D_Kernels_tmp = (G1 + G2 - 2*G)/2
    D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    D_Kernels[[i]] = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
  }

  return(D_Kernels)

}

# compute the single kernel
"dist2.cimlr" = function( x, c = NA ) {
  # set the parameters for x
  if(is.na(c)) {
    c = x
  }

  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }

  # compute the distance
  dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) +
    (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) -
    2 * (x%*%t(c))

  dist[which(dist<0,arr.ind=TRUE)] = 0

  return(dist)

}

# normalizes a symmetric kernel
"dn.cimlr" = function( w, type ) {

  # compute the sum of any column
  w = w * dim(w)[1]
  D = apply(abs(w),MARGIN=1,FUN=sum)

  # type "ave" returns D^-1*W
  if(type=="ave") {
    D = 1 / D
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% w
  }
  else {
    stop("Invalid type!")
  }

  return(wn)

}

# umkl function
"umkl.cimlr2" = function( D, beta = NA, u = 100) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  print(u)
  u = u
  logU = log(u)

  # compute Hbeta
  res_hbeta = Hbeta(D, beta)
  H = res_hbeta$H
  thisP = res_hbeta$P

  betamin = -Inf
  betamax = Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin = beta
      if(abs(betamax)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamax) / 2
      }
    }
    else {
      betamax = beta
      if(abs(betamin)==Inf) {
        beta = beta / 2
      }
      else {
        beta = (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }

  return(thisP)

}

# umkl function
"umkl.cimlr.gl" = function( D, beta = NA, u = 50 ) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = u
  logU = log(u)

  # compute Hbeta
  res_hbeta = Hbeta(D, beta)
  H = res_hbeta$H
  thisP = res_hbeta$P

  betamin = -Inf
  betamax = Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin = beta
      if(abs(betamax)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamax) / 2
      }
    }
    else {
      betamax = beta
      if(abs(betamin)==Inf) {
        beta = beta / 2
      }
      else {
        beta = (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }

  return(list(thisP = thisP, beta = beta))

}




# todo:----
# 1. 把cimlr的function移过来；本地也source 相同的scipt
# 2. 改变cimlr的函数, 除了umkl以外其他参数可以外给
# 3. 给完全相同参数，看pkd和cimlr的差别
# fitting functions ----
multi_kernel_weight = function(cur_sim){
  sim_data_weight = list(data_type1 = list(data = cur_sim$data1, weight = cur_sim$weight1),
                         data_type2 = list(data = cur_sim$data2, weight = cur_sim$weight2),
                         data_type3 = list(data = cur_sim$data3, weight = cur_sim$weight3))
  unlist(lapply(sim_data_weight,function(data_type){
    multi_kernels_gl(x = data_type$data, feature_weight = data_type$weight)
  }))
}

propose_kernel_dist_fit = function(mk_list,truelabel = NA, cur_sim= NA, par1, is_tsne = F,c = 4){
  res = dfl_simlr(kernel_list = mk_list,
                  c=c,
                  alpha0 = par1$alpha0,
                  alpha = par1$alpha,
                  rho = par1$rho,
                  beta = par1$beta,
                  gamma = par1$gamma,
                  normalize = 0,
                  print_details = F,
                  umkl = F,
                  is_tsne = is_tsne)
  if(is.na(truelabel)) {truelabel = cur_sim$truelabel}
  cluster = spectralClustering(as.matrix(res$S), K = c)
  nmi_tmp=igraph::compare(truelabel, cluster, method = "nmi")
  rand_tmp=igraph::compare(truelabel, cluster, method = "rand")
  silhouette_tmp =
    summary(silhouette_similarity(cluster,as.matrix(res$S)))$avg.width

  return(list(res = res,
              evaluation = c(nmi = nmi_tmp,
                             rand = rand_tmp,
                             silhouette = silhouette_tmp,
                             method = "dfl_simlr",
                             w1 = sum(res$w[1:22]),
                             w2 = sum(res$w[23:44]),
                             w3 = sum(res$w[45:66])
              )))
}

propose_kernel_fit = function(mk_list,cur_sim= NA, truelabel = NA,par1, c = 4){
  res = fl_simlr(kernel_list = mk_list,
                 c=c,
                 alpha0 = par1$alpha0,
                 alpha = par1$alpha,
                 rho = par1$rho,
                 beta = par1$beta,
                 gamma = par1$gamma,
                 normalize = 0,
                 print_details = F,
                 umkl = F)
  if(is.na(truelabel)) {truelabel = cur_sim$truelabel}
  cluster = spectralClustering(as.matrix(res$S), K = c)
  nmi_tmp=igraph::compare(truelabel, cluster, method = "nmi")
  rand_tmp=igraph::compare(truelabel, cluster, method = "rand")
  silhouette_tmp =
    summary(silhouette_similarity(cluster,as.matrix(res$S)))$avg.width

  return(list(res = res,
              evaluation = c(nmi = nmi_tmp,
                             rand = rand_tmp,
                             silhouette = silhouette_tmp,
                             method = "fl_simlr_w",
                             w1 = sum(res$w[1:22]),
                             w2 = sum(res$w[23:44]),
                             w3 = sum(res$w[45:66])
              )))
}

absnf_fit = function(cur_sim= NA, truelabel = NA,c = 4){
  sim_data_weight = list(data_type1 = list(data = cur_sim$data1, weight = cur_sim$weight1),
                         data_type2 = list(data = cur_sim$data2, weight = cur_sim$weight2),
                         data_type3 = list(data = cur_sim$data3, weight = cur_sim$weight3))
  aff_mat_list = lapply(sim_data_weight, function(data_type) affinityMatrix(dist2_w(as.matrix(data_type$data),
                                                                                    as.matrix(data_type$data),
                                                                                    weight = data_type$weight)))
  S = SNF(aff_mat_list, K = 20)
  cluster = spectralClustering(S, K = c)
  if(is.na(truelabel)) {truelabel = cur_sim$truelabel}
  nmi_tmp=igraph::compare(truelabel, cluster, method = "nmi")
  rand_tmp=igraph::compare(truelabel, cluster, method = "rand")
  silhouette_tmp =
    summary(silhouette_similarity(cluster,S))$avg.width
  return(list(res = list(S = S, cluster = cluster),
              evaluation = c(nmi = nmi_tmp,
                             rand = rand_tmp,
                             silhouette = silhouette_tmp,
                             method = "absnf",
                             w1 = NA,
                             w2 = NA,
                             w3 = NA
              )))
}

snf_fit = function(cur_sim= NA, truelabel = NA,c = 4){
  aff_mat_list = lapply(cur_sim[1:3], function(y) affinityMatrix(dist2(as.matrix(y),as.matrix(y))))
  S = SNF(aff_mat_list, K = 20)
  cluster = spectralClustering(S, K = c)
  if(is.na(truelabel)) {truelabel = cur_sim$truelabel}
  nmi_tmp=igraph::compare(truelabel, cluster, method = "nmi")
  rand_tmp=igraph::compare(truelabel, cluster, method = "rand")
  silhouette_tmp =
    summary(silhouette_similarity(cluster,S))$avg.width
  return(list(res = list(S = S, cluster = cluster),
              evaluation = c(nmi = nmi_tmp,
                             rand = rand_tmp,
                             silhouette = silhouette_tmp,
                             method = "snf",
                             w1 = NA,
                             w2 = NA,
                             w3 = NA
              )))
}

cimlr_fit = function(cur_sim= NA, truelabel = NA,c = 4){
  res = CIMLR(lapply(cur_sim[1:3], function(y) t(y)), c = c)
  if(is.na(truelabel)) {truelabel = cur_sim$truelabel}
  nmi_tmp=igraph::compare(truelabel, res$y$cluster, method = "nmi")
  rand_tmp=igraph::compare(truelabel, res$y$cluster, method = "rand")
  silhouette_tmp =
    summary(silhouette_similarity(res$y$cluster,res$S))$avg.width
  return(list(res = res, evaluation = c(nmi = nmi_tmp,
                                        rand = rand_tmp,
                                        silhouette = silhouette_tmp,
                                        method = "cimlr",
                                        w1 = sum(res$alphaK[1:55]),
                                        w2 = sum(res$alphaK[56:110]),
                                        w3 = sum(res$alphaK[111:165])
  )))
}

# further check ----
multi_kernels_gl_p = function (x, sigma =seq(2, 1, -0.1), cores.ratio = 0.5,p = 2, standardize = 1) {
  # data _normalization
  # x = apply(x, 2, function(y) {
  #   (y-min(y))/(max(y)-min(y))
  # })

  N = dim(x)[1]
  KK = 0
  mu = seq(0.1,1, 0.1)
  allk = c(10,30)

  ## calculate distance
  if(standardization == 1){
    x = standardNormalization(x)
  }else if(standardization == 2){
    x = apply(x, 2, function(y) {(y-min(y))/(max(y)-min(y))})}


  # distance calc and sort

  Diff = as.matrix(dist(x_fun,method = "minkowski", p = p))


  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    library(Matrix)
    library(SNFtool)
    D_kernels = function(Diff_fun = NA,
                         allk_fun = 30,
                         sigma_fun = 2,
                         is.normal = T,
                         is.SP = T,
                         standardization = standardize) {
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
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, x_fun = x, Diff_fun = Diff, B = allk,
                                                                  sig_df_list = sigma, mu_aff_list = mu, KK_fun = KK) {
    if (B[l] < (nrow(x_fun) - 1)) {
      for (j in 1:length(sig_df_list)){
        name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
        D_Kernels[[name]] = Matrix(D_kernels(Diff_fun = Diff_fun,
                                             allk_fun = B[l],
                                             sigma_fun = sig_df_list[[j]],
                                             standardization = 0))}
    }



    return(D_Kernels)
  }))

  stopCluster(cl)
  return(D_Kernels)
}


###### SEPARATE STEP FUNCTIONS##### ----
# calc mk ----

# input:
# * data_list: m data types with n_samples*p_feature structure
# * alpha: parameter of local term

# Notes:
# 1. if normalized_type is not 0, then we shouldn't do the data normalization'
# 2. Data normalization is outsiede of the kernel calc function

Dist_kernel_generation = function(sim,
                                  alpha = 0,
                                  k = 30,
                                  sigma = 2,
                                  allk = 30,
                                  is.square = T,
                                  data_used = 1:3,
                                  dist_used = NA,
                                  local_term = F,
                                  diffusion_term = F,
                                  mutual_knn = F,
                                  normalized_type = "no", # options: no_coef, no, gph, ave
                                  trans_type = "P",
                                  t=2,
                                  kernel_f = NA
                                  ){
  # calc kernel ----

  if(is.na(kernel_f)){ ## old kernel generation
    with_coef = ifelse(normalized_type == "no_coef", F, T)
    if("dist_list" %in% names(sim)){
      Kernel = unlist(lapply(sim$dist_list, function(x) multi_kernels_gl(Diff = x,
                                                                         is.square = is.square,
                                                                         with_coef = with_coef,
                                                                         sigma =sigma,
                                                                         allk = allk)))
      ## dist kerenl only look at single kernel right now
    }else{
      if(length(sim)<6){ data_used = 1}
      Kernel = unlist(lapply(sim[data_used], function(x) multi_kernels_gl(x = as.matrix(x),
                                                                          is.square = is.square,
                                                                          with_coef = with_coef,
                                                                          sigma = sigma,
                                                                          allk = allk)))
    }
  }else{ # new kernel taking different forms
    if("dist_list" %in% names(sim)){
      Kernel = unlist(lapply(sim$dist_list, function(x) multi_kernels_2(Diff = x,
                                                                        is.square = is.square,
                                                                        sigma =sigma,
                                                                        allk = allk,
                                                                        kernel_f = kernel_f)))
      ## dist kerenl only look at single kernel right now
    }else{
      if(length(sim)<6){ data_used = 1}
      Kernel = unlist(lapply(sim[data_used], function(x) multi_kernels_2(x = as.matrix(x),
                                                                         is.square = is.square,
                                                                         sigma = sigma,
                                                                         allk = allk,
                                                                         kernel_f = kernel_f)))
    }
  }
  # local enhance ----
  if(alpha>0 & diffusion_term){
    Kernel_diffusion = lapply(Kernel, function(K){
      n_samp = nrow(K)
      if(mutual_knn){
        K_knn = mutual_knn_func(s = K,k = k)
        K_knn = (K_knn+t(K_knn))/2
      }else{
        K_knn_init = dominateset(K,k+1)
        K_knn = (K_knn_init+t(K_knn_init))/2
      }
      # diag(K_knn) = 0
      if(trans_type == "P"){
        transition_mat = dn.cimlr(K_knn, type="ave")
      }else if (trans_type == "T"){
        transition_mat = transition.fields(K_knn)
      }

      if(is.numeric(t[1])){
        L_list = list()
        for(j in 1:length(t)){
          if(t[j] == 0){
            L = K_knn
          }else{
            transition_t = diag(1, nrow(K_knn))
            for(i in 1:t[j]){
              transition_t = transition_t %*% t(transition_mat)
            }
            L = K %*% transition_t
          }
          L_list[[j]] = L
        }
        L = Reduce("+", L_list)/length(t)
        L = (L+t(L))/2
      }else if (t == "sym"){
        L = transition_mat %*% K %*% transition_mat
      }else if(t == "infinity"){
        eps = 2e-16
        eig = eigen(transition_mat)
        U <- eig$vectors
        d <- eig$values - eps
        if(trans_type == "T"){
          d <- (1 - alpha) * d/(1 - alpha * d^2)
          D <- diag(d)
          L <- U %*% D %*% t(U)
        }else{
          d = 1/(1-(1-alpha)*d)
          D <- diag(d)
          L_trans <- U %*% D %*% t(U)
          L = alpha* K %*% L_trans
        }

      }

      if(t!="infinity") {L = L * alpha}

      L
    })

    if(t!="infinity"){
      Kernel = lapply(Kernel, function(x) (1-alpha)*x)
      Kernel = Map("+",Kernel,Kernel_diffusion)
    }else{
      Kernel = Kernel_diffusion
    }

  }

  # normalizetion ----
  if(normalized_type %in% c("ave","gph")) { # kernel normaliztion
    Kernel = lapply(Kernel, function(x){x = dn(x, normalized_type)})
  }

  # kernel_dist generation
  Dk = lapply(Kernel,function(x){ # kernel to dist kernel
    dist_kernels(as.matrix(x))
  })

  # dk_tibble output ----
  if("dist_list" %in% names(sim)){
    dk_tibble = tibble(distance = map_chr(str_split(names(Dk),pattern  = "\\."),1),
                       dist_kernel = Dk,
                       kernel = Kernel
                       ) %>%
      group_by(distance) %>%
      nest() %>%
      mutate(dist_kernel = map(data, "dist_kernel"),
             kernel = map(data,"kernel")) %>%
      dplyr::select(-c(data)) %>%
      ungroup()
    dk_tibble = rbind(dk_tibble, tibble(distance = "all",
                                  dist_kernel = list(Dk),
                                  kernel = list(Kernel)))
    # dk_tibble = tibble(distance = c("BC","HM","Weighted","Unweigted","all"),
    #                    dist_kernel = list(Dk$BC,
    #                                       Dk$HM,
    #                                       Dk$Weighted,
    #                                       Dk$Unweighted,
    #                                       Dk
    #                    ),
    #                    kernel = list(Kernel$BC,
    #                                  Kernel$HM,
    #                                  Kernel$Weighted,
    #                                  Kernel$Unweighted,
    #                                  Kernel
    #                    ))
  }else{
    dk_tibble = tibble(distance = "Euclidean",
                       dist_kernel = list(Dk),
                       kernel = list(Kernel))
  }


  return(dk_tibble)
}
# calc initialization ----
initialization = function(D_Kernels,
                          c = 4,# number of clusters
                          k = 30,# #neighbors
                          cores.ratio = 1,
                          network_diffusion = T){

  if(!is.null(dim(D_Kernels))) D_Kernels = list(D_Kernels)
  num = ncol(D_Kernels[[1]]) # samples
  r = -1
  beta = 0.8
  ## weight init ----
  alphaK = 1/rep(length(D_Kernels), length(D_Kernels))
  names(alphaK) = names(D_Kernels)
  print(k)
  ## S init ----
  distX = array(0, c(dim(D_Kernels[[1]])[1], dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX/length(D_Kernels)
  S0 = max(max(distX)) - distX
  if(network_diffusion){S0 = network.diffusion(S0, k)}
  S0 = as.matrix(S0)
  S0 = dn.cimlr(S0, "ave")
  S = (S0 + t(S0))/2
  # sort dist
  res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x,
                                                               index.return = TRUE)))
  distX1 = array(0, c(nrow(distX), ncol(distX)))
  idx = array(0, c(nrow(distX), ncol(distX)))
  for (i in 1:nrow(distX)) {
    distX1[i, ] = res[[i]]$x
    idx[i, ] = res[[i]]$ix
  }
  ## H init ----
  D0 = diag(apply(S, MARGIN = 2, FUN = sum))
  L0 = D0 - S
  eig1_res = eig1(L0, c, 0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  F_eig1 = dn.cimlr(F_eig1, "ave")
  ## par estimation init ----
  A = array(0, c(num, num))
  di = distX1[, 2:(k + 2)]
  rr = 0.5 * (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1,
                                      FUN = sum))
  id = idx[, 2:(k + 2)]
  numerator = (apply(array(0, c(length(di[, k + 1]), dim(di)[2])),
                     MARGIN = 2, FUN = function(x) {
                       x = di[, k + 1]
                     }) - di)
  temp = (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1, FUN = sum) +
            .Machine$double.eps)
  denominator = apply(array(0, c(length(temp), dim(di)[2])),
                      MARGIN = 2, FUN = function(x) {
                        x = temp
                      })
  temp = numerator/denominator
  a = apply(array(0, c(length(t(1:num)), dim(di)[2])), MARGIN = 2,
            FUN = function(x) {
              x = 1:num
            })
  A[cbind(as.vector(a), as.vector(id))] = as.vector(temp)
  if (r <= 0) {
    r = mean(rr)
  }
  lambda = max(mean(rr), 0)
  A[is.nan(A)] = 0

  ## return res ----
  cat(lambda,"\n", r, "\n")
  print(S0[1:5,1:5])
  return(list(alphaK = alphaK, # weight init
              F_eig1 = F_eig1, # H init
              evs_eig1 = evs_eig1,
              S = S, # S init
              S0 = S0,
              distX = distX,
              distX1 = distX1, # distX sorted by row
              idx = idx, # distX sort index of each row
              lambda = lambda,
              r = r
  ))
}
# iteration ----

### !!!!!output EACH iteration result!!!!
CIMLR_optimization = function(D_Kernels,
                              Kernel = NA,
                              alphaK, # weight init
                              F_eig1, # H init
                              S, # S init
                              S0,
                              evs_eig1,
                              distX,
                              distX1, # distX sorted by row
                              idx, # distX sort index of each row
                              lambda,
                              r,
                              c = 4,# number of clusters
                              k = 30,# #neighbors
                              cores.ratio = 1,
                              alpha = 0,
                              update_neighb = F,
                              u = 150,
                              network_diffusion = F){
  if(!is.null(dim(D_Kernels))) D_Kernels = list(D_Kernels)
  if(!is.null(dim(Kernel))) Kernel = list(Kernel)
  NITER = 30 # max iteration
  num = ncol(D_Kernels[[1]]) # samples
  beta = 0.8
  no.dim = c
  # start iteration ----
  S_iter = list()
  converge = vector()
  for(iter in 1:NITER){
    ## S ----
    cat("Iteration: ", iter, "\n")
    distf = L2_distance_1(t(F_eig1), t(F_eig1))
    A = array(0, c(num, num))
    b = idx[, 2:dim(idx)[2]]
    a = apply(array(0, c(num, ncol(b))), MARGIN = 2, FUN = function(x) {
      x = 1:num
    })
    inda = cbind(as.vector(a), as.vector(b))
    ad = (distX[inda] + lambda * distf[inda])/2/r
    dim(ad) = c(num, ncol(b))
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx", c_input, c_output))
    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    S = (1 - beta) * A + beta * S

    if(network_diffusion){
      S = network.diffusion(S, k)
    }

    S = as.matrix(S)
    S = (S + t(S))/2
    S_iter[[iter]] = S
    ### update local
    if(alpha>0 & !is.na(Kernel) & update_neighb){
      neighb_ind = neighbor_index_generate(S,k)
      D_Kernels = update_Dk(Kernel,neighb_ind,alpha)
    }
    ## H ----
    D = diag(apply(S, MARGIN = 2, FUN = sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L, c, 0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    F_eig1 = dn.cimlr(F_eig1, "ave")
    F_eig1 = (1 - beta) * F_old + beta * F_eig1
    # evs_eig1 = cbind(evs_eig1, ev_eig1)
    ## w ----
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps + D_Kernels[[i]]) *
        (S + .Machine$double.eps)
      DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
    }
    alphaK0 = umkl.cimlr2(DD,u = u)
    alphaK0 = alphaK0/sum(alphaK0)
    alphaK = (1 - beta) * alphaK + beta * alphaK0
    alphaK = alphaK/sum(alphaK)
    print(alphaK)
    names(alphaK) = names(D_Kernels)
    ## converge judge and par update ----
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c + 1)])
    converge[iter] = fn2 - fn1
    if (iter < 10) {
      if (ev_eig1[length(ev_eig1)] > 1e-06) {
        lambda = 1.5 * lambda
        r = r/1.01
      }
    }
    else {
      if (converge[iter] > 1.01 * converge[iter - 1]) {
        S = S_old
        if (converge[iter - 1] > 0.2) {
          warning("Maybe you should set a larger value of c.")
        }
        break
      }
    }
    ## update weighted sum ----
    S_old = S
    distX = D_Kernels[[1]] * alphaK[1]
    if(length(D_Kernels)>1){
      for (i in 2:length(D_Kernels)) {
        distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
      }
    }

    res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x,
                                                                 index.return = TRUE)))
    distX1 = array(0, c(nrow(distX), ncol(distX)))
    idx = array(0, c(nrow(distX), ncol(distX)))
    for (i in 1:nrow(distX)) {
      distX1[i, ] = res[[i]]$x
      idx[i, ] = res[[i]]$ix
    }
  }

  # output ----
  LF = F_eig1
  D = diag(apply(S, MARGIN = 2, FUN = sum))
  L = D - S
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  if (length(no.dim) == 1) {
    U_index = seq(ncol(U), (ncol(U) - no.dim + 1))
    F_last = U[, U_index]
  }

  cat("Performing Kmeans.\n")
  y = kmeans(F_last, c, nstart = 200)
  ydata = tsne(S)
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["S_iter"]] = S_iter
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["converge"]] = converge
  results[["LF"]] = LF
  return(results)
}




solQP_optimization = function(D_Kernels = D_Kernels,
                              Kernel = NA,
                              alphaK = alphaK, # weight init
                              F_eig1 = F_eig1, # H init
                              S = S, # S init
                              S0 = S0,
                              distX = distX,
                              distX1 = distX1, # distX sorted by row
                              idx = idx, # distX sort index of each row
                              lambda = lambda,
                              r = r,
                              c = 4,# number of clusters
                              k = 30,# #neighbors
                              cores.ratio = 1,
                              umkl = T,
                              new_beta = F,
                              u = 150,
                              alpha = 0,
                              update_neighb = F,
                              print_details = F,
                              stopping = 1e-5,
                              diag_sum = NULL,
                              network_diffusion = F){
  old_w = alphaK
  old_S = S
  old_L = F_eig1
  beta = 0.8
  rho = 1 # using umkl in the optimization, rho here only for the optimization value calc
  n_sam = ncol(S)
  beta = lambda
  gamma = r
  opt_value_vec = NULL
  n_ite = 30
  S_iter = list()
  # Start iteration ----
  for (iter in 1:n_ite) {
    cat(iter)

    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = - distX + gamma * (old_L %*% t(old_L))
    I_n = diag(1, n_sam)

    ### Update S ----
    for (i in 1:n_sam) {
      constraint = constraint_produce(i, diag_sum,S0, n_sam, I_n)
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = constraint$Amat,
                            bvec = constraint$bvec,
                            meq = constraint$meq)
      new_S[,i] = QP_results$solution
    }

    if(network_diffusion){
      diag(new_S) = 0
      new_S = (1-beta)*new_S+beta*old_S
      new_S = network.diffusion(new_S,k)
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric
    S_iter[[iter]] = new_S
    ### update local
    if(alpha>0 & !is.na(Kernel) & update_neighb){
      neighb_ind = neighbor_index_generate(S,k)
      D_Kernels = update_Dk(Kernel,neighb_ind,alpha)
    }
    # if(network.diffusion){new_S = network.diffusion(new_S, K = k)}

    ### Update L ----
    Laplacian = diag(apply(new_S,1, sum)) - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w ----
    first_term = vector() # 1st term in optimization
    n_ker = length(D_Kernels)
    for (ll in 1:n_ker) {
      first_term[ll] = sum(-D_Kernels[[ll]] * new_S)
    }
    if(umkl){
      DD = vector()
      for (i in 1:length(D_Kernels)) {
        temp = (.Machine$double.eps + D_Kernels[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD,u = u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      for (ll in 1:n_ker) {
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    distX = matrix(0,n_sam,n_sam)
    for(i in 1:length(new_w)){
      distX = distX + new_w[i]*D_Kernels[[i]]
    }



    ### calculate optimization value ----
    opt_value = sum(distX * new_S) +
      # rho * sum(new_w * log(new_w)) +
      beta * norm(new_S, type = 'F')^2 +
      gamma * sum(diag(t(new_L) %*% (diag(apply(new_S, 1,sum)) - new_S) %*% new_L))
    opt_value_vec = c(opt_value_vec, opt_value)


    ### Whether to stop ----
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }
  cluster = kmeans(new_L, c, nstart = 200)$cluster

  res = list(w = new_w,
             S = new_S,
             S_iter = S_iter,
             L = new_L,
             cluster = cluster,
             beta_umkl = beta_umkl,
             opt_value_vec = opt_value_vec)
  return(res)
}

solQP_kernel_optimization = function(Kernel = NA,
                              alphaK = alphaK, # weight init
                              F_eig1 = F_eig1, # H init
                              S = S, # S init
                              S0 = S0,
                              distX = distX,
                              distX1 = distX1, # distX sorted by row
                              idx = idx, # distX sort index of each row
                              lambda = lambda,
                              r = r,
                              c = 4,# number of clusters
                              k = 30,# #neighbors
                              cores.ratio = 1,
                              umkl = T,
                              new_beta = F,
                              u = 150,
                              alpha = 0,
                              update_neighb = F,
                              print_details = F,
                              stopping = 1e-5,
                              diag_sum = NULL,
                              network_diffusion = F){
  old_w = alphaK
  old_S = S
  old_L = F_eig1
  beta = 0.8
  rho = 1 # using umkl in the optimization, rho here only for the optimization value calc
  n_sam = ncol(S)
  beta = lambda
  gamma = r
  opt_value_vec = NULL
  n_ite = 30
  S_iter = list()
  weighted_sum_kernels = matrix(0, n_sam, n_sam)
  for(l in 1:length(Kernel)){
    weighted_sum_kernels =
      weighted_sum_kernels + as.matrix(Kernel[[l]])*old_w[[l]]
  }
  # Start iteration ----
  for (iter in 1:n_ite) {
    cat(iter)

    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = weighted_sum_kernels + gamma * (old_L %*% t(old_L))
    I_n = diag(1, n_sam)

    ### Update S ----
    for (i in 1:n_sam) {
      constraint = constraint_produce(i, diag_sum,S0, n_sam, I_n)
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = constraint$Amat,
                            bvec = constraint$bvec,
                            meq = constraint$meq)
      new_S[,i] = QP_results$solution
    }

    if(network_diffusion){
      diag(new_S) = 0
      new_S = (1-beta)*new_S+beta*old_S
      new_S = network.diffusion(new_S,k)
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric
    S_iter[[iter]] = new_S
    ### update local
    if(alpha>0 & !is.na(Kernel) & update_neighb){
      neighb_ind = neighbor_index_generate(S,k)
      D_Kernels = update_Dk(Kernel,neighb_ind,alpha)
    }

    ### Update L ----
    Laplacian = diag(apply(new_S,1, sum)) - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w ----
    first_term = vector() # 1st term in optimization
    n_ker = length(Kernel)
    for (ll in 1:n_ker) {
      first_term[ll] = sum(Kernel[[ll]] * new_S)
    }
    if(umkl){
      DD = vector()
      for (i in 1:length(Kernel)) {
        temp = (.Machine$double.eps + -Kernel[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD,u = u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      for (ll in 1:n_ker) {
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }

    weighted_sum_kernels = matrix(0, n_sam, n_sam)
    for(l in 1:length(Kernel)){
      weighted_sum_kernels =
        weighted_sum_kernels + as.matrix(Kernel[[l]])*old_w[[l]]
    }



    ### calculate optimization value ----
    opt_value = sum(-weighted_sum_kernels * new_S) +
      # rho * sum(new_w * log(new_w)) +
      beta * norm(new_S, type = 'F')^2 +
      gamma * sum(diag(t(new_L) %*% (diag(apply(new_S, 1,sum)) - new_S) %*% new_L))
    opt_value_vec = c(opt_value_vec, opt_value)


    ### Whether to stop ----
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }
  cluster = kmeans(new_L, c, nstart = 200)$cluster

  res = list(w = new_w,
             S = new_S,
             S_iter = S_iter,
             L = new_L,
             cluster = cluster,
             beta_umkl = beta_umkl,
             opt_value_vec = opt_value_vec)
  return(res)
}



# fitting method ----

propose_kernel_dist_cimlr_fit = function(sim, # currently taking 2 kinds, results from simulation_3 and microb_simulation
                                         c = 4, # [FIXME]: number of clusters, will change to 3 automaticlly if the sim is from microb_simulation
                                         k = 30, # number of neighbours
                                         network_diffusion = F, # original SIMLR method, whether to perform network diffusion for the learned S in each step
                                         # parameter for kernel generation
                                         data_used = 1:3, # determine which data in sim to integrate
                                         dist_used = "single_and_all", # used for simulation obj from microb_simulation, if single and all, look at results from single kernel and integ all dist, otherwise take all dist only
                                         alpha = 0, # parameter for term L
                                         diffusion_term = F, # whether to include diffustion term L with alpha
                                         is.square = F, # if TRUE, follow SIMLR algorithm with d^4 as nominator
                                         mutual_knn = F, # if TRUE, the affinity matrix entry = 1 if i,j are mutual neighbours
                                         normalized_type = "no", # normalization type for the affinity mat, "ave" means row-wise normalization, "gph" means row-and-column-wise normalization
                                         t = 1, # determine how many steps will the kernel diffusion go. If NA, perform MKM,otw, perform K(M^T)^t
                                         trans_type = "P", # if P, D^{-1}A, if T, follow def in Wang B(2018)
                                         local_term = F, # [no use now],indicate whether to use local term,
                                         update_neighb = F, # [no use now], parameter for local term with updated neighbors,
                                         sigma = 2,
                                         allk = 30,
                                         u = 150,
                                         c_single = NA,
                                         kernel_f = NA
                                         ){
  dk_tibble = Dist_kernel_generation(sim,
                                     data_used = data_used,
                                     dist_used = dist_used,
                                     alpha = alpha,
                                     k = k,
                                     normalized_type = normalized_type,
                                     trans_type = trans_type,
                                     mutual_knn = mutual_knn,
                                     is.square = is.square,
                                     local_term = local_term,
                                     diffusion_term = diffusion_term,
                                     t = t,
                                     sigma = sigma,
                                     allk = allk,
                                     kernel_f = kernel_f
                                     )
  if(dist_used == "all"){
    dk_tibble = dk_tibble %>%
      filter(distance %in% c("all","Euclidean"))
  }

  res_tibble = dk_tibble %>%
    mutate(res = pmap(list(dk = dist_kernel,
                           kernel = kernel,
                           dist = distance),
                      function(dk,kernel, dist){
                        print(names(dk))
                        print(c)
                        print(dist)
                        if(dist != "all" & dist!= "Euclidean"){ c = c_single}
                        print(c)
                        init = initialization(dk, c = c, k = k,network_diffusion = network_diffusion)
                        res = CIMLR_optimization(D_Kernels = dk,
                                                 Kernel = kernel,
                                                 alphaK = init$alphaK, # weight init
                                                 F_eig1 = init$F_eig1, # H init
                                                 S = init$S, # S init
                                                 S0 = init$S0,
                                                 distX = init$distX,
                                                 distX1 = init$distX1, # distX sorted by row
                                                 idx = init$idx, # distX sort index of each row
                                                 lambda = init$lambda,
                                                 r = init$r,
                                                 c = c,# number of clusters
                                                 k = k,# #neighbors
                                                 cores.ratio = 1,
                                                 alpha = alpha,
                                                 update_neighb = update_neighb,
                                                 u = u,
                                                 network_diffusion = network_diffusion)
                        return(res)

    })) %>%
    dplyr::select(-c(dist_kernel,kernel))
  return(res_tibble)
}


propose_kernel_dist_sQP_fit = function(sim,data_used = 1:3,
                                       alpha = 0,
                                       diag_sum = "initial",
                                       c,
                                       update_neighb = F,
                                       network_diffusion = F){
  if(update_neighb){
    kernel_list = Dist_kernel_generation(sim[data_used])
  }else{
    kernel_list = Dist_kernel_generation(sim[data_used],alpha = alpha)
  }
  init = initialization(kernel_list$Dk)
  res = solQP_optimization(D_Kernels = kernel_list$Dk,
                           Kernel = kernel_list$Kernel,
                           alphaK = init$alphaK, # weight init
                           F_eig1 = init$F_eig1, # H init
                           S = init$S, # S init
                           S0 = init$S0,
                           distX = init$distX,
                           distX1 = init$distX1, # distX sorted by row
                           idx = init$idx, # distX sort index of each row
                           lambda = init$lambda,
                           r = init$r,
                           c = 4,# number of clusters
                           k = 30,# #neighbors
                           cores.ratio = 1,
                           alpha = alpha,
                           update_neighb = update_neighb,
                           diag_sum = diag_sum,
                           u = length(data_used)*50,
                           network_diffusion = network_diffusion)
  return(res)
}

propose_kernel_sQP_fit = function(sim,data_used = 1:3,
                                       alpha = 0,
                                       diag_sum = "initial",
                                       c,
                                       update_neighb = F,
                                       network_diffusion = F){
  if(update_neighb){
    kernel_list = Dist_kernel_generation(sim[data_used])
  }else{
    kernel_list = Dist_kernel_generation(sim[data_used],alpha = alpha)
  }
  init = initialization(kernel_list$Kernel)
  res = solQP_kernel_optimization(
                           Kernel = kernel_list$Kernel,
                           alphaK = init$alphaK, # weight init
                           F_eig1 = init$F_eig1, # H init
                           S = init$S, # S init
                           S0 = init$S0,
                           distX = init$distX,
                           distX1 = init$distX1, # distX sorted by row
                           idx = init$idx, # distX sort index of each row
                           lambda = init$lambda,
                           r = init$r,
                           c = 4,# number of clusters
                           k = 30,# #neighbors
                           cores.ratio = 1,
                           alpha = alpha,
                           update_neighb = update_neighb,
                           diag_sum = diag_sum,
                           u = length(data_used)*50,
                           network_diffusion = network_diffusion)
  return(res)
}

SNF_fitting = function(sim,
                       c = 4,
                       k = 30,
                       c_def = T,
                       kernel_f = "SNF",
                       NUMC = 2:11){
  # construct affinity mat
  if(kernel_f == "SNF"){
    if("dist_list" %in% names(sim)){
      if(c_def){c = 3}
      aff_list = lapply(sim$dist_list, affinityMatrix)
      distance = "all"
    }else{
      if("data3" %in% names(sim)){
        dist_list = lapply(sim[1:3], dist2)
      }else{
        dist_list = list(dist2(as.matrix(sim$data)))
      }

      aff_list = lapply(dist_list, affinityMatrix)

    }
  }else{
    dk_tib = Dist_kernel_generation(sim = sim,
                                    sigma = 2,
                                    allk = 30,
                                    kernel_f = kernel_f)
    if("dist_list" %in% names(sim)){
      aff_list = lapply(unlist(dk_tib$kernel[1:4]), as.matrix)
    }else{
      aff_list = dk_tib$kernel[[1]]
    }
  }
  distance = "Euclidean"
  # fusion
  if(length(aff_list)>1){
    W = SNF(aff_list)
  }else(
    W = aff_list[[1]]
  )
  if(is.na(c)){
    c = estimateNumberOfClustersGivenGraph(W, NUMC = NUMC)[[1]]
  }
  clusters = spectralClustering(W, c)

  if(!is.null(sim$truelabel)){
    nmi = compare(clusters,sim$truelabel,"nmi")
  }else{
    nmi = 0
  }

  # calculate cluster, construct result
  res= list()
  res[["cluster"]] = clusters
  res[["S"]] = W
  res_tibble = tibble(distance = distance,
                      res = list(res),
                      nmi = nmi
                      #adjust_rand = igraph::compare(spectralClustering(W, c), sim$truelabel,"adjusted.rand"),
                      )
  return(res_tibble)
}

## microbiome fitting ----

# this function directly using the distance kernel as input for CIMLR method
single_all_kernel_distance_cimlr_compare = function(data,
                                                    network_diffusion = T,
                                                    alpha = 0,
                                                    update_neighb = F){

  # fitting with single and multiple distance kernel ----
  compare_tib = tibble(distance = names(data$dist_list),
                       distance_data = data$dist_list,
                       truth = list(data$truelabel)) # create distance kernel compare tibble

  compare_tib = rbind(compare_tib,
                      tibble(distance = "all",
                             distance_data = list(data$dist_list),
                             truth = list(data$truelabel))) %>% # add kernel and kernel_distance
    mutate(nest_kernel = map2(distance,distance_data, function(dist_type,dist){
      if(dist_type == "all"){
        kernel_gernerated = lapply(dist, function(x) D_kernels(Diff_fun = x))
        kernel_distance = unlist(lapply(kernel_gernerated, function(x) x$D_Kernels))
        kernel = unlist(lapply(kernel_gernerated, function(x) x$Kernels))
        nest_kernel = tibble(kernel_distance = list(kernel_distance),
                             kernel = list(kernel))
      }else{
        kernel_gernerated = D_kernels(Diff_fun = dist)
        kernel_distance = kernel_gernerated$D_Kernels
        kernel = kernel_gernerated$Kernels
        nest_kernel = tibble(kernel_distance = list(kernel_distance),
                             kernel = list(kernel))
      }
    })) %>%
    unnest(nest_kernel) %>% # perform initialization
    mutate(init = map(kernel_distance, function(dk){
      initialization(D_Kernels = dk,c = 3,k = 30 )
    })) %>% # Perform optimization using SIMLR
    mutate(res_cim = map2(init, kernel_distance, function(init, dk){
      res_cim = CIMLR_optimization(D_Kernels = dk,
                                   alphaK = init$alphaK, # weight init
                                   F_eig1 = init$F_eig1, # H init
                                   S = init$S, # S init
                                   S0 = init$S0,
                                   distX = init$distX,
                                   distX1 = init$distX1, # distX sorted by row
                                   idx = init$idx, # distX sort index of each row
                                   lambda = init$lambda,
                                   r = init$r,
                                   c = 3,# number of clusters
                                   k = 30,# #neighbors，
                                   u = 3,
                                   alpha = alpha,
                                   update_neighb = update_neighb,
                                   network_diffusion = network_diffusion,
                                   cores.ratio = 1)
      return(res_cim)
    }))

  compare_tib = compare_tib %>%
    mutate(cluster = map(res_cim, function(x) x$y$cluster),
           weight = map(res_cim, function(x) x$alphaK)) %>%
    dplyr::select(-c(res_cim, distance_data, kernel_distance, kernel, init)) %>%
    mutate(method = "kernel dist cimlr")


  compare_tib = compare_tib %>%
    mutate(nmi = map2_dbl(cluster,truth, function(x,y){compare(x,y,method = "nmi")}),
           adjust_rand = map2_dbl(cluster,truth, function(x,y){compare(x,y,method = "adjusted.rand")})
    )

  return(compare_tib)

}

# utils ----

dn = function( w, type ) {

  # compute the sum of any column
  D = apply(w,MARGIN=2,FUN=sum)

  # type "ave" returns D^-1*W
  if(type=="ave") {
    D = 1 / D
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% w
  }
  # type "gph" returns D^-1/2*W*D^-1/2
  else if(type=="gph") {
    D = 1 / sqrt(D)
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% (w %*% D)
  }
  else {
    stop("Invalid type!")
  }

  return(wn)

} # used to normalize the similarity matrix

mutual_knn_func = function(s = NA,k){
  n = nrow(s) # sample size
  knn = dominateset(s,KK = k+1)
  diag(knn) = 0
  mut_knn = sqrt(knn*t(knn))
  # then enhance weak connectivity
  # first identify points with less than k/2 neighb
  isolated_id = (1:n)[apply(mut_knn,1,function(x) sum(x!=0)<(k/2))]
  # create connections between these points with: 1. first k/2 neighb 2. the last k/2 neighb with connection < k/2
  for(i in isolated_id){
    sort_i = sort(s[i,], decreasing = T, index.return = T)$ix
    # find points to connect with i: low connection and first k
    neighb_and_low_conn_id = intersect(sort_i[1:k], isolated_id)
    # connect
    mut_knn[i,neighb_and_low_conn_id] = mut_knn[neighb_and_low_conn_id,i] = s[i, neighb_and_low_conn_id]
    # connect the shortest path to avoid 0 connetion
    mut_knn[i,sort_i[2]] = s[i, sort_i[2]]
    }
  mut_knn
}

constraint_produce = function(i,
                              diag_sum_sign,
                              S0 = NA,
                              n_sam,
                              I_n){
  if(is.null(diag_sum_sign)){
    A_i = t(rbind(rep(1,n_sam), I_n))
    b = c(1, rep(0,n_sam))
    meq = 1
  }else{
    diag_sign = rep(0,n_sam)
    diag_sign[i] = 1
    A_i = t(rbind(rep(1,n_sam), diag_sign,I_n))
    meq = 2
    if(diag_sum_sign == "initial"){
      b = c(1, diag(S0)[i],rep(0,n_sam))
    }else if(diag_sum_sign == 0){
      b = c(1,0,rep(0,n_sam))
    }
  }
  return(list(Amat = A_i,
              bvec = b,
              meq = meq))
}
neighbor_index_generate = function(S, k = 30){
  sorting = apply(S, MARGIN = 1, function(x) sort(x, index.return = T,decreasing = T))
  idx = list()
  for(i in 1:length(sorting)){
    idx[[i]] = sorting[[i]]$ix[2:(k+1)]
  }
  idx
}

update_Dk =
  function(Kernel, neighbor_ind, alpha = 0){
    if(alpha>0){
      Kernel_local = lapply(Kernel, function(x){
        x_l = matrix(0, nrow(x), ncol(x))
        for(i in 1:nrow(x)){
          x_l[i, neighbor_ind[[i]]] = x[i, neighbor_ind[[i]]]
        }
        x_l*alpha
      } )
      Kernel = Map("+",Kernel,Kernel_local)
    }
    Dk = lapply(Kernel,dist_kernels)
    Dk
  }

"umkl.cimlr2" = function( D, beta = NA, u = 100) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = u
  logU = log(u)

  # compute Hbeta
  res_hbeta = Hbeta(D, beta)
  H = res_hbeta$H
  thisP = res_hbeta$P

  betamin = -Inf
  betamax = Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin = beta
      if(abs(betamax)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamax) / 2
      }
    }
    else {
      betamax = beta
      if(abs(betamin)==Inf) {
        beta = beta / 2
      }
      else {
        beta = (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }

  return(thisP)

}




