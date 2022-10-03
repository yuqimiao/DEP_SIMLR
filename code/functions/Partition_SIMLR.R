partition_SIMLR = function(sim,
                           data_used = 1:3,
                           k=20,
                           rho = 0.001,
                           c = 4,
                           c_single = 2,
                           c_iterative = F,
                           numc = 2:11,
                           network_diffusion = F,
                           diff_term = F,
                           alpha = 0,
                           kernel_normalize = "no_coef",
                           Y_normalize = F,
                           weight_type =  "consensus",
                           sigma = 2,
                           allk = 30,
                           kernel_f = NA
                           ){
  # estimate number of clusters ----
  if(is.na(c_single)){
    c_single = map_dbl(sim[data_used], function(x){
      score = Est_num_clus(X = x, NUMC = numc)$K1
      numc[which.min(score)]
    })
  }

  if(is.na(c)){
    if(!c_iterative){
      score = Est_num_clus(sim[data_used], NUMC = numc)$K1
      c = numc[which.min(score)]
    }
  }

  # calc distance kernels ----
  dk_tib = Dist_kernel_generation(sim = sim,
                                  data_used = data_used,
                                  diffusion_term = diff_term,
                                  alpha = alpha,
                                  sigma = sigma,
                                  allk = allk,
                                  normalized_type = kernel_normalize,
                                  kernel_f = kernel_f)


  if(nrow(dk_tib) == 1){
    distX_ls = dk_tib$dist_kernel[[1]]
  }else{
    dk_tib = dk_tib %>% filter(distance == "all")
    distX_ls = dk_tib$dist_kernel[[1]]
  }
  # initialization ----
  S = length(distX_ls)
  if(length(c_single)==1){c_single = rep(c_single,S)}
  n = nrow(distX_ls[[1]])

  w0 = rep(1/S,S)
  w_cur = w0



  r = -1
  # eta = 0.8 # learning rate
  # lambda_set = 1

  ## Get Z0
  Z0 = lapply(distX_ls, function(x){
    x = as.matrix(x)
    max(x)-x
  })
  Z_cur = Z0
  if(network_diffusion){Z0 = lapply(Z0, function(x) network.diffusion(x, K = k))}
  initial_list = lapply(distX_ls, function(distX){
    print(dim(distX))
    # rank every row of the kernel distance
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x # ith row is the ranked ith row of kernel distance
      idx[i,] = res[[i]]$ix # index of the ranks of ith row
    }
    A = array(0,c(n,n))
    di = distX1[,2:(k+2)]
    rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum)) # parameter estimation
    id = idx[,2:(k+2)]

    if(r<=0) {
      r = mean(rr)
    }
    lambda = max(mean(rr),0)
    # lambda = lambda_set

    return(list(distX = distX,
                distX1=distX1,
                idx = idx,
                lambda = lambda))
  })
  ## Get F0
  cat("THINK CLUSTER NUMBER ESTIMATION FOR SINGLE DATA")

  F0 = lapply(1:length(Z0), function(i){
    L = diag(rowSums(Z0[[i]]))-Z0[[i]]
    F_eig1 = eig1(L, c = c_single[i], isMax = 0)$eigvec
    # F_eig1 = dn.cimlr(F_eig1, "ave")
    F_eig1
  })

  weight_L0 = matrix(0,n,n)
  for(s in 1:S){
    L = diag(1,n)-F0[[s]]%*%t(F0[[s]])*2
    weight_L0 = weight_L0+w0[[s]]*L
  }
  Y0 = eig1(weight_L0, isMax = 0, c =c )$eigvec
  Z_cur = Z0
  F_cur = F0
  w_cur=w0
  Y_cur = Y0

  # interation ----
  converge = 100
  converge_cur = converge
  for(t in 1:30){
    cat("iteration", t, "\n")
    # update Z ----
    Z_pre = Z_cur
    for(s in 1:S){
      # updata each data type separately
      F_eig1 = F_cur[[s]]
      distX = initial_list[[s]]$distX
      distX1 = initial_list[[s]]$distX1
      idx = initial_list[[s]]$idx
      lambda = initial_list[[s]]$lambda
      r = initial_list[[s]]$lambda
      distf = L2_distance_1(t(F_eig1),t(F_eig1))
      A = array(0,c(n,n))
      b = idx[,2:dim(idx)[2]]
      a = apply(array(0,c(n,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:n })
      inda = cbind(as.vector(a),as.vector(b)) # rank of each row aligned
      ad = (distX[inda]+lambda*distf[inda])/2/r
      dim(ad) = c(n,ncol(b))

      # call the c function for the optimization
      c_input = -t(ad)
      c_output = t(ad)
      ad = t(.Call("projsplx_R",c_input,c_output))
      A[inda] = as.vector(ad)
      A[is.nan(A)] = 0
      A = (A + t(A)) / 2
      # Z_cur[[s]] = (1 - eta) * Z_cur[[s]] + eta * A
      Z_cur[[s]] = A
      if(network_diffusion){Z_cur[[s]] = network.diffusion(Z_cur[[s]], k)}
      # Z_cur[[s]] = dn(Z_cur[[s]],"ave")
    }
    # update F ----
    F_pre = F_cur
    for(s in 1:S){
      # updata each data type separately
      L = diag(rowSums(Z_cur[[s]]))-Z_cur[[s]]
      lambda = initial_list[[s]]$lambda
      r = initial_list[[s]]$lambda
      F_eig1 = eig1(lambda*L+rho*w_cur[[s]]*(diag(1, n)-2*Y_cur%*%t(Y_cur)), isMax = 0, c =c_single[[s]] )$eigvec
      # F_eig1 = dn.cimlr(F_eig1, "ave")
      # F_eig1 = (1 - eta) * F_pre[[s]] + eta * F_eig1
      F_cur[[s]] = F_eig1

    }

    trace_single_cur = unlist(lapply(1:S, function(s){
      sum(diag(t(F_cur[[s]])%*% (diag(apply(Z_cur[[s]],1,sum))-Z_cur[[s]]) %*%F_cur[[s]]))
    }))
    # update Y----
    # if(t>1) {Y_pre = Y_cur}
    w_pre = w_cur
    weight_L = matrix(0,n,n)
    for(s in 1:S){
      L = diag(1,n)-F_cur[[s]]%*%t(F_cur[[s]])*2
      weight_L = weight_L+w_cur[[s]]*L
    }
    if(is.na(c)){
      X_int = list.cbind(map2(w_cur, F_cur, function(w, Fc){
        w*Fc
      }))
      kernel_dist_X_int = as.matrix(multi_kernels_2(x = X_int,
                                          allk = 30,
                                          sigma = 2)[[1]])
      score = Est_num_clus(distance = kernel_dist_X_int, NUMC = numc)$K1
      c = numc[which.min(score)]
    }
    Y_cur = eig1(weight_L, isMax = 0, c = c )$eigvec
    # Y_cur = dn.cimlr(Y_cur,"ave")
    # update weight ----
    w_pre = w_cur
    if(weight_type == "consensus"){
      YYt = Y_cur %*% t(Y_cur)
      for(s in 1:S){
        FsFst = F_cur[[s]]%*%t(F_cur[[s]])
        w_cur[[s]] = 1/sqrt(2*sum(abs(YYt-FsFst)))
      }

      w_cur = w_cur/sum(w_cur)
      print(w_cur)
    }else if (weight_type == "trace"){
      w_cur = (1/trace_single_cur)/sum((1/trace_single_cur))
    }

    # calculate optimization function ----
    converge_pre = converge_cur
    converge_cur = 0
    for(s in 1:S){
      distX = initial_list[[s]]$distX
      beta = gamma = initial_list[[s]]$lambda
      term1 = sum(distX * Z_cur[[s]])
      term2 = sum(Z_cur[[s]]^2) * beta
      term3 = sum(diag( t(F_cur[[s]])%*% (diag(apply(Z_cur[[s]],1,sum))-Z_cur[[s]]) %*%F_cur[[s]])) * gamma
      term4 = rho* w_cur[[s]] * sum((Y_cur %*% t(Y_cur)-F_cur[[s]]%*%t(F_cur[[s]]))^2)
      converge_cur = converge_cur+term1+term2+term3+term4
      # cat(term1,term2,term3,term4, "\n")
    }
    converge = abs(converge_cur-converge_pre)
    print(converge)
    if(converge<1e-10){
      break
    }

  }

  if(Y_normalize){Y_cur = dn.cimlr(Y_cur,"ave")}
  # get results ----
  res_S = Y_cur %*% t(Y_cur)
  # clusters = spectralClustering(res_S, K = 4)
  clusters = kmeans(Y_cur,c = c, nstart = 200)$cluster
  if(!is.null(sim$truelabel)){
    nmi = compare(clusters,sim$truelabel,"nmi")
  }else{
    nmi = 0
  }


  return(list(Y_cur = Y_cur,
              F_cur = F_cur,
              c_single = c_single,
              trace_single = trace_single_cur,
              Z_cur = Z_cur,
              clusters = clusters,
              nmi = nmi,
              w_cur = w_cur))

}


partition_SIMLR_large = function(distX_ls,
                                 k=20,
                                 rho = 0.001,
                                 c = 3,
                                 c_single =3,
                                 truelabel
){

  S = length(distX_ls)
  print(S)
  n = nrow(distX_ls[[1]])

  w0 = rep(1/S,S)
  w_cur = w0


  if(length(c_single)==1){c_single = rep(c_single,S)}
  r = -1
  # eta = 0.8 # learning rate
  # lambda_set = 1

  ## Get Z0
  Z0 = lapply(distX_ls, function(x){
    x = as.matrix(x)
    max(x)-x
  })
  Z_cur = Z0
  initial_list = lapply(distX_ls, function(distX){
    print(dim(distX))
    # rank every row of the kernel distance
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x # ith row is the ranked ith row of kernel distance
      idx[i,] = res[[i]]$ix # index of the ranks of ith row
    }
    A = array(0,c(n,n))
    di = distX1[,2:(k+2)]
    rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum)) # parameter estimation
    id = idx[,2:(k+2)]

    if(r<=0) {
      r = mean(rr)
    }
    lambda = max(mean(rr),0)
    # lambda = lambda_set

    return(list(distX = distX,
                distX1=distX1,
                idx = idx,
                lambda = lambda))
  })
  ## Get F0
  cat("THINK CLUSTER NUMBER ESTIMATION FOR SINGLE DATA")

  F0 = lapply(1:length(Z0), function(i){
    L = diag(rowSums(Z0[[i]]))-Z0[[i]]
    F_eig1 = eig1(L, c = c_single[i], isMax = 0)$eigvec
    # F_eig1 = dn.cimlr(F_eig1, "ave")
    F_eig1
  })

  weight_L0 = matrix(0,n,n)
  for(s in 1:S){
    L = diag(1,n)-F0[[s]]%*%t(F0[[s]])*2
    weight_L0 = weight_L0+w0[[s]]*L
  }
  Y0 = eig1(weight_L0, isMax = 0, c =c )$eigvec
  Z_cur = Z0
  F_cur = F0
  w_cur=w0
  Y_cur = Y0

  converge = 100
  converge_cur = converge
  for(t in 1:30){
    cat("iteration", t, "\n")
    # update Z ----
    Z_pre = Z_cur
    for(s in 1:S){
      # updata each data type separately
      F_eig1 = F_cur[[s]]
      distX = initial_list[[s]]$distX
      distX1 = initial_list[[s]]$distX1
      idx = initial_list[[s]]$idx
      lambda = initial_list[[s]]$lambda
      r = initial_list[[s]]$lambda
      distf = L2_distance_1(t(F_eig1),t(F_eig1))
      A = array(0,c(n,n))
      b = idx[,2:dim(idx)[2]]
      a = apply(array(0,c(n,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:n })
      inda = cbind(as.vector(a),as.vector(b)) # rank of each row aligned
      ad = (distX[inda]+lambda*distf[inda])/2/r
      dim(ad) = c(n,ncol(b))

      # call the c function for the optimization
      c_input = -t(ad)
      c_output = t(ad)
      ad = t(.Call("projsplx_R",c_input,c_output))
      A[inda] = as.vector(ad)
      A[is.nan(A)] = 0
      A = (A + t(A)) / 2
      # Z_cur[[s]] = (1 - eta) * Z_cur[[s]] + eta * A
      Z_cur[[s]] = A
      # Z_cur[[s]] = network.diffusion(Z_cur[[s]], k)
      # Z_cur[[s]] = dn(Z_cur[[s]],"ave")
    }
    # update F ----
    F_pre = F_cur
    for(s in 1:S){
      # updata each data type separately
      L = diag(rowSums(Z_cur[[s]]))-Z_cur[[s]]
      lambda = initial_list[[s]]$lambda
      r = initial_list[[s]]$lambda
      F_eig1 = eig1(lambda*L+rho*w_cur[[s]]*(diag(1, n)-2*Y_cur%*%t(Y_cur)), isMax = 0, c =c_single[[s]] )$eigvec
      # F_eig1 = dn.cimlr(F_eig1, "ave")
      # F_eig1 = (1 - eta) * F_pre[[s]] + eta * F_eig1
      F_cur[[s]] = F_eig1

    }
    # update Y----
    # if(t>1) {Y_pre = Y_cur}
    w_pre = w_cur
    weight_L = matrix(0,n,n)
    for(s in 1:S){
      L = diag(rowSums(F_cur[[s]]%*%t(F_cur[[s]])))-F_cur[[s]]%*%t(F_cur[[s]])*2
      weight_L = weight_L+w_cur[[s]]*L
    }
    Y_cur = eig1(weight_L, isMax = 0, c =c )$eigvec
    # Y_cur = dn.cimlr(Y_cur,"ave")
    # update weight ----
    w_pre = w_cur
    YYt = Y_cur %*% t(Y_cur)
    for(s in 1:S){
      FsFst = F_cur[[s]]%*%t(F_cur[[s]])
      w_cur[[s]] = 1/sqrt(2*sum(abs(YYt-FsFst)))
    }

    w_cur = w_cur/sum(w_cur)
    print(w_cur)

    # calculate optimization function ----
    converge_pre = converge_cur
    converge_cur = 0
    for(s in 1:S){
      distX = initial_list[[s]]$distX
      beta = gamma = initial_list[[s]]$lambda
      term1 = sum(distX * Z_cur[[s]])
      term2 = sum(Z_cur[[s]]^2) * beta
      term3 = sum(diag( t(F_cur[[s]])%*% (diag(apply(Z_cur[[s]],1,sum))-Z_cur[[s]]) %*%F_cur[[s]])) * gamma
      term4 = rho* w_cur[[s]] * sum((Y_cur %*% t(Y_cur)-F_cur[[s]]%*%t(F_cur[[s]]))^2)
      converge_cur = converge_cur+term1+term2+term3+term4
      # cat(term1,term2,term3,term4, "\n")
    }
    converge = abs(converge_cur-converge_pre)
    print(converge)
    if(converge<1e-10){
      break
    }

  }

  # Y_cur = dn.cimlr(Y_cur,"ave")
  # res_S = Z_cur
  # clusters = spectralClustering(res_S, K = 4)
  clusters = kmeans(Y_cur,c = c, nstart = 200)$cluster
  nmi = compare(clusters,truelabel,"nmi")

  return(list(Y_cur = Y_cur,
              Z_cur = Z_cur,
              clusters = clusters,
              nmi = nmi,
              w_cur = w_cur))

}


# utils ----
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
