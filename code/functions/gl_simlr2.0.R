# Header ----
library(parallel)
library(clValid)
library(dplyr)
library(SNFtool)
library(igraph)
library(Matrix)
library(quadprog)
## glsimlr( no adaptive learning) ----
gl_simlr = function(data_list = NA, dist_list = NA, kernel_list = NA, B = 30, c = 4, max = 30, tol = 1e-10,
                    rho = 0.1,gamma = 0.5,beta = 0.5,delta = 0,kernel_type = "Dk",sigma_dk = 2,
                    mu_aff = 0.5, isdk = F, standardization = 2){

  # kernel calculation
  if(is.na(kernel_list)){
    if(is.na(dist_list)){
      if(kernel_type == "Dk"){
        dist_kernel_list = lapply(data_list, FUN = function(x) D_kernels(x_fun = x,allk_fun = B,sigma_fun = sigma_dk,standardization = standardization) )
        if(isdk){
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})
        }else{
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
        }

      }else if(kernel_type == "affinity"){
        if(standardization == 1){
          data_list = lapply(data_list, standardNormalization)
        }else if(standardization == 2){
          data_list = lapply(data_list, function(x) {
            apply(x, 2, function(y) {
              (y-min(y))/(max(y)-min(y))
            })
          })
        }
        kernel_list = lapply(data_list, FUN = function(x) affinityMatrix(dist2(x)^(1/2),sigma = mu_aff))
      }
    }else{ ## didn't test whether applicable to dist kernel
      if(kernel_type == "Dk"){
        dist_kernel_list = lapply(data_list, FUN = function(x) D_kernels(Diff_fun = x,allk_fun = B,sigma_fun = sigma_dk,standardization = standardization))
        if(isdk){
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})
        }else{
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
        }
      }else if(kernel_type == "affinity"){
        kernel_list = lapply(dist_list, FUN = function(x) affinityMatrix(as.matrix(x)^(1/2)))
      }
    }
  }

  # kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
  local_kernel_list = lapply(kernel_list, FUN = function(x) dominateset(x, KK = B))
  # D_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})

  # Kernel normalization
  norm_kernel_list = lapply(kernel_list, FUN = function(x){
    D = diag(rowSums(x))
    solve(D^(1/2))%*%as.matrix(x)%*%solve(D^(1/2))
  }
  )

  norm_local_kernel_list = lapply(local_kernel_list, FUN = function(x){
    D = diag(rowSums(x))
    solve(D^(1/2))%*%as.matrix(x)%*%solve(D^(1/2))
  }
  )
  # Global local term calculation
  global_term = kernel_similarity(norm_kernel_list)
  local_term = kernel_similarity(norm_local_kernel_list)

  # w0, H0 initialization
  M = length(kernel_list)
  n = dim(kernel_list[[1]])[1]
  w = rep(1/M,M)
  K = matrix(0,n,n)


  for(i in 1:M){
    K = K+norm_kernel_list[[i]]*w[i]
  }

  D = diag(rowSums(K))
  L = D-K
  L = solve(D^(1/2))%*%L%*%solve(D^(1/2))
  eig_res = eig1(L, c, isMax = 0)
  H = eig_res$eigvec

  # iteration
  criteria = NULL
  w_all = NULL
  trace_term = NULL
  for(t in 1:max){
    ## iteration of new w given H
    trace = NULL
    for(i in 1:M){
      L_i = diag(rowSums(norm_kernel_list[[i]])) - norm_kernel_list[[i]]
      trace = c(trace, sum(diag(t(H)%*%L_i%*%H)))
      w[i]=exp((delta*(global_term[i]+beta*local_term[i])-gamma*sum(diag(t(H)%*%L_i%*%H)))/rho)
    }
    w = w/sum(w)
    w_all = rbind(w_all, w)
    trace_term = rbind(trace_term, trace)

    ## iteration of H given w
    K = matrix(0,n,n)
    for(i in 1:M){
      K = K+norm_kernel_list[[i]]*w[i]
    }
    D = diag(rowSums(K))
    L = D-K
    L = solve(D^(1/2))%*%L%*%solve(D^(1/2))
    eig_res = eig1(L, c, isMax = 0)
    H = eig_res$eigvec

    ## stop  critieria
    obj = -delta*(sum(w*global_term+w*beta*local_term))+gamma*sum(diag(t(H)%*%L%*%H))+rho*sum(w*log(w))
    criteria = c(criteria,obj)
    if(t > 10){
      if ((criteria[t]-criteria[t-1])<tol){
        cluster = kmeans(H,c,nstart = 200)$cluster
        break
      }
    }

  }
  return(list(w = w,
              K = K,
              w_all = w_all,
              trace_term = trace_term,
              global_term = global_term,
              local_term = local_term,
              cluster = cluster,
              iteration = t))
}

# fl simlr
fl_simlr = function(data_list = NA, # a list of data matrices to integrate
                     kernel_list = NA,
                     k = 30, # #neighbors
                     c = 3, # #clusters
                     kernel_fun = "D_kernels", # the definition of kernels
                     rho = 0.1, # tuning parameter for entropy penalty
                     alpha0,# tuning parameter for <S, K>
                     alpha, # tuning parameter for <S, W>
                     beta = NA, # tuning parameter for ||S||
                     gamma = NA, # tuning parameter for laplacian
                     stopping = 10^-3, # stopping rule
                     n_ite = 50, # max number of iterations
                     normalize = 2,
                     print_details = F, # whether to print stepwise details
                     self_weight = F,
                     umkl = F,
                     u = 50,
                     center_kernel = F
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }

  if(center_kernel){
    s_Kmat_ls = lapply(s_Kmat_ls, center_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, center_kernel)
  }

  if(umkl){
    Dk_ls = lapply(s_Kmat_ls, dist_kernels)
    Dq_ls = lapply(s_Qmat_ls, dist_kernels)
  }

  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    w_s_Kmat = w_s_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
      w_s_Qmat = w_s_Qmat + old_w[ll] * s_Qmat_ls[[ll]]
    }

    # parameter estimation
    if(is.na(beta)&is.na(gamma)){
      D_w_Kmat = dist_kernels(w_s_Kmat)
      D_w_Qmat = dist_kernels(w_s_Qmat)
      par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
        s_x = sort(x)
        rep(s_x[k+2],k) - s_x[2:(k+1)]
      }))
      print(par_est)
      gamma = beta = par_est
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * w_s_Kmat + alpha * w_s_Qmat + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric

    ### Update L
    Laplacian = I_n - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
      DD = vector()
      for (i in 1:length(Dk_ls)) {
        temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD, u=u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*s_Kmat_ls[[ll]]+alpha*s_Qmat_ls[[ll]]) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  beta_umkl = ifelse(umkl, umkl_res$beta, "NULL")
  return(list(w = new_w,
              S = new_S,
              L = new_L,
              beta_umkl = beta_umkl))
}


# dfl simlr ----
dfl_simlr = function(data_list = NA, # a list of data matrices to integrate
                     kernel_list = NA,
                     k = 30, # #neighbors
                     c = 3, # #clusters
                     kernel_fun = "D_kernels", # the definition of kernels
                     rho = 0.1, # tuning parameter for entropy penalty
                     alpha0,# tuning parameter for <S, K>
                     alpha, # tuning parameter for <S, W>
                     beta = NA, # tuning parameter for ||S||
                     gamma = NA, # tuning parameter for laplacian
                     stopping = 10^-3, # stopping rule
                     n_ite = 50, # max number of iterations
                     normalize = 2,
                     print_details = F, # whether to print stepwise details
                     self_weight = F,
                     umkl = F,
                     u = 50,
                     center_kernel = F,
                     new_beta = F,
                     network.diffusion = F,
                     is_tsne = F # further cluster using tsne
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }

  if(center_kernel){
    s_Kmat_ls = lapply(s_Kmat_ls, center_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, center_kernel)
  }

  Dk_ls = lapply(s_Kmat_ls, dist_kernels)
  Dq_ls = lapply(s_Qmat_ls, dist_kernels)


  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  if(network.diffusion){old_S = network.diffusion(old_S, K = k)}
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    w_s_Kmat = w_s_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
      w_s_Qmat = w_s_Qmat + old_w[ll] * s_Qmat_ls[[ll]]
    }
    # weighted distance kernel
    D_w_Kmat = dist_kernels(w_s_Kmat)
    D_w_Qmat = dist_kernels(w_s_Qmat)

    # parameter estimation
    if(is.na(beta)){
      if(new_beta){
        D_Lmat = dist2(old_L)
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat+gamma*D_Lmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }else{
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }
      print(par_est)
      beta = par_est
    }

    if(is.na(gamma)){
      gamma = beta
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * (-D_w_Kmat) + alpha * (-D_w_Qmat) + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric
    if(network.diffusion){new_S = network.diffusion(new_S, K = k)}

    ### Update L
    Laplacian = I_n - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
      DD = vector()
      for (i in 1:length(Dk_ls)) {
        temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD,u = u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*-(Dk_ls[[ll]])+alpha*-(Dk_ls[[ll]])) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  beta_umkl = ifelse(umkl, umkl_res$beta,"NULL")

  if(is_tsne){
    S = as.matrix(new_S)
    new_L = tsne(X = S, k = c, initial_config = new_L)
  }
  cluster = kmeans(new_L, c, nstart = 200)$cluster
  return(list(w = new_w,
              S = new_S,
              L = new_L,
              cluster = cluster,
              beta_umkl = beta_umkl))
}

dfl_simlr_mod = function(data_list = NA, # a list of data matrices to integrate
                     kernel_list = NA,
                     k = 30, # #neighbors
                     c = 3, # #clusters
                     kernel_fun = "D_kernels", # the definition of kernels
                     rho = 0.1, # tuning parameter for entropy penalty
                     alpha0,# tuning parameter for <S, K>
                     alpha, # tuning parameter for <S, W>
                     beta = NA, # tuning parameter for ||S||
                     gamma = NA, # tuning parameter for laplacian
                     stopping = 10^-3, # stopping rule
                     n_ite = 50, # max number of iterations
                     normalize = 2,
                     print_details = F, # whether to print stepwise details
                     self_weight = F,
                     umkl = F,
                     u = 50,
                     center_kernel = F,
                     new_beta = F,
                     network.diffusion = F,
                     is_tsne = F # further cluster using tsne
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }

  if(center_kernel){
    s_Kmat_ls = lapply(s_Kmat_ls, center_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, center_kernel)
  }

  Dk_ls = lapply(s_Kmat_ls, dist_kernels)
  Dq_ls = lapply(s_Qmat_ls, dist_kernels)


  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  if(network.diffusion){old_S = network.diffusion(old_S, K = k)}
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(diag(apply(old_S,1,sum))-old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    D_w_Kmat =D_w_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      D_w_Kmat = D_w_Kmat + old_w[ll] * Dk_ls[[ll]]
      D_w_Qmat = D_w_Qmat + old_w[ll] * Dq_ls[[ll]]
    }
    # weighted distance kernel
    # D_w_Kmat = dist_kernels(w_s_Kmat)
    # D_w_Qmat = dist_kernels(w_s_Qmat)

    # parameter estimation
    if(is.na(beta)){
      if(new_beta){
        D_Lmat = dist2(old_L)
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat+gamma*D_Lmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }else{
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }
      print(par_est)
      beta = par_est
    }

    if(is.na(gamma)){
      gamma = beta
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * (-D_w_Kmat) + alpha * (-D_w_Qmat) + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric
    if(network.diffusion){new_S = network.diffusion(new_S, K = k)}

    ### Update L
    Laplacian = diag(apply(new_S,1,sum)) - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
      DD = vector()
      for (i in 1:length(Dk_ls)) {
        temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD,u = u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*-(Dk_ls[[ll]])+alpha*-(Dq_ls[[ll]])) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  beta_umkl = ifelse(umkl, umkl_res$beta,"NULL")

  if(is_tsne){
    S = as.matrix(new_S)
    new_L = tsne(X = S, k = c, initial_config = new_L)
  }
  cluster = kmeans(new_L, c, nstart = 200)$cluster
  return(list(w = new_w,
              S = new_S,
              L = new_L,
              cluster = cluster,
              beta_umkl = beta_umkl))
}

dfl_simlr_Sneighbor = function(data_list = NA, # a list of data matrices to integrate
                               kernel_list = NA,
                               k = 30, # #neighbors
                               c = 3, # #clusters
                               kernel_fun = "D_kernels", # the definition of kernels
                               rho = 0.1, # tuning parameter for entropy penalty
                               alpha0,# tuning parameter for <S, K>
                               alpha, # tuning parameter for <S, W>
                               beta = NA, # tuning parameter for ||S||
                               gamma = NA, # tuning parameter for laplacian
                               stopping = 10^-3, # stopping rule
                               n_ite = 50, # max number of iterations
                               normalize = 2,
                               print_details = F, # whether to print stepwise details
                               self_weight = F,
                               umkl = F,
                               u = 50,
                               center_kernel = F,
                               new_beta = F,
                               network.diffusion = F,
                               is_tsne = F # further cluster using tsne
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  # s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      # s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    # s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }

  if(center_kernel){
    s_Kmat_ls = lapply(s_Kmat_ls, center_kernel)
    # s_Qmat_ls = lapply(s_Qmat_ls, center_kernel)
  }

  Dk_ls = lapply(s_Kmat_ls, dist_kernels)
  # Dq_ls = lapply(s_Qmat_ls, dist_kernels)


  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  if(network.diffusion){old_S = network.diffusion(old_S, K = k)}
  neigbor_index = apply(old_S, 1, function(x) sort(x,index.return = T))
  Dq_ls = lapply(Dk_ls, function(x){

  })
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(diag(apply(old_S,1,sum))-old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    D_w_Kmat =D_w_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      D_w_Kmat = D_w_Kmat + old_w[ll] * Dk_ls[[ll]]
      D_w_Qmat = D_w_Qmat + old_w[ll] * Dq_ls[[ll]]
    }
    # weighted distance kernel
    # D_w_Kmat = dist_kernels(w_s_Kmat)
    # D_w_Qmat = dist_kernels(w_s_Qmat)

    # parameter estimation
    if(is.na(beta)){
      if(new_beta){
        D_Lmat = dist2(old_L)
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat+gamma*D_Lmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }else{
        par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
          s_x = sort(x)
          rep(s_x[k+2],k) - s_x[2:(k+1)]
        }))
      }
      print(par_est)
      beta = par_est
    }

    if(is.na(gamma)){
      gamma = beta
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * (-D_w_Kmat) + alpha * (-D_w_Qmat) + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric
    if(network.diffusion){new_S = network.diffusion(new_S, K = k)}

    ### Update L
    Laplacian = diag(apply(new_S,1,sum)) - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
      DD = vector()
      for (i in 1:length(Dk_ls)) {
        temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      umkl_res = umkl.cimlr.gl(DD,u = u)
      new_w = umkl_res$thisP
      beta_umkl = umkl_res$beta
      new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*-(Dk_ls[[ll]])+alpha*-(Dq_ls[[ll]])) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  beta_umkl = ifelse(umkl, umkl_res$beta,"NULL")

  if(is_tsne){
    S = as.matrix(new_S)
    new_L = tsne(X = S, k = c, initial_config = new_L)
  }
  cluster = kmeans(new_L, c, nstart = 200)$cluster
  return(list(w = new_w,
              S = new_S,
              L = new_L,
              cluster = cluster,
              beta_umkl = beta_umkl))
}





# tsne ----
"tsne" <- function(X, initial_config = NULL, k = 2, max_iter = 1000, min_cost = 0, epoch = 100) {

  cat("Performing t-SNE.\n")

  momentum = 0.8
  final_momentum = 0.8
  mom_switch_iter = 250

  epsilon = 500
  min_gain = 0.01
  initial_P_gain = 4

  n = nrow(X)

  eps = .Machine$double.eps

  if (!is.null(initial_config) && is.matrix(initial_config)) {
    print(dim(initial_config))
    print(n)
    print(k)
    if (nrow(initial_config) != n | ncol(initial_config) != k) {
      stop('initial_config argument does not match necessary configuration for X')
    }
    ydata = initial_config
    initial_P_gain = 1

  }
  else {
    ydata = matrix(rnorm(k * n),n)
  }

  P = X
  P = 0.5 * (P + t(P))

  P[P < eps]<-eps
  P = P/sum(P)

  P = P * initial_P_gain
  grads = matrix(0,nrow(ydata),ncol(ydata))
  incs = matrix(0,nrow(ydata),ncol(ydata))
  gains = matrix(1,nrow(ydata),ncol(ydata))

  for (iter in 1:max_iter) {

    if (iter %% epoch == 0) {
      cost = sum(apply(P * log((P+eps)/(Q+eps)),1,sum))

      cat("Epoch: Iteration #",iter," error is: ",cost,"\n")

      if (cost < min_cost) {
        break
      }
    }

    sum_ydata = apply((ydata^2),1,sum)
    num =  1/(1 + sum_ydata + sweep(-2*ydata %*% t(ydata),2,-t(sum_ydata)))
    diag(num) = 0
    Q = num / sum(num)

    if (any(is.nan(num))) {
      message ('NaN in grad. descent')
    }

    Q[Q < eps] = eps

    stiffnesses = (P-Q) * num
    grads = 4 * (diag(apply(stiffnesses,2,sum)) - stiffnesses) %*% ydata

    gains = (gains + .2) * abs(sign(grads) != sign(incs)) + gains * .8 * abs(sign(grads) == sign(incs))
    gains[gains < min_gain] = min_gain

    incs = momentum * incs - epsilon * (gains * grads)
    ydata = ydata + incs
    ydata = sweep(ydata,2,apply(ydata,2,mean))

    # we are constraining the ydata
    ydata[ydata < -100] = -100
    ydata[ydata > 100] = 100

    if (iter == mom_switch_iter) {
      momentum = final_momentum
    }

    if (iter == 100 && is.null(initial_config)) {
      P = P/4
    }

  }

  return(ydata)

}
