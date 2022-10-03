source("code/R/network.diffusion.R")
source("code/functions/Partition_CIMLR_2.0.R")
CIMLR_kernel = function (kernel_list, c, no.dim = NA, k = 10, cores.ratio = 1)
{
  if (is.na(no.dim)) {
    no.dim = c
  }

  NITER = 30
  num = ncol(kernel_list[[1]])
  r = -1
  beta = 0.8
  D_Kernels = lapply(kernel_list, function(kernel) dist_kernels(kernel))
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
    ad = t(.Call("projsplx_R", c_input, c_output))
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
      temp = (.Machine$double.eps + D_Kernels[[i]]) *
        (S + .Machine$double.eps)
      DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
    }
    alphaK0 = umkl.cimlr(DD)
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

  cluster = spectralClustering(S, K=c)

  results = list(cluster = cluster,
                 S = S)

  return(results)
}

# umkl function
"umkl.cimlr" = function( D, beta = NA ) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = 150
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
"Hbeta" = function( D, beta ) {

  D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
  P = exp(-D * beta)
  sumP = sum(P)
  H = log(sumP) + beta * sum(D * P) / sumP
  P = P / sumP

  return(list(H=H,P=P))

}
