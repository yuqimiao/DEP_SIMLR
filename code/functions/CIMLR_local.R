# library(SIMLR)
library(Matrix)
library("parallel")
library(gplots)
library(igraph)
library("clValid")

source("./code/R/SIMLR.R")
source("./code/R/compute.multiple.kernel.cimlr.R")
source("./code/R/network.diffusion.R")
source("./code/R/utils.simlr.R")
source("./code/R/tsne.R")

dyn.load("./code/R/projsplx_R.so")

## removing techs: tsne, network.diffusion
CIMLR_rmTechs = function (X, c, no.dim = NA, k = 10, cores.ratio = 1)
{
  if (is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  NITER = 30
  num = ncol(X[[1]])
  r = -1
  beta = 0.8
  cat("Computing the multiple Kernels.\n")
  for (data_types in 1:length(X)) {
    curr_X = X[[data_types]]
    if (data_types == 1) {
      D_Kernels = multiple.kernel.cimlr(t(curr_X), cores.ratio)
    }
    else {
      D_Kernels = c(D_Kernels, multiple.kernel.cimlr(t(curr_X),
                                                     cores.ratio))
    }
  }
  # initialization ----
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
  # S0 = network.diffusion(S0, k)
  S0 = as.matrix(S0)
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

  # iteration ----
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
    # S = network.diffusion(S, k)
    S = as.matrix(S)
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
  if (length(no.dim) == 1) {
    U_index = seq(ncol(U), (ncol(U) - no.dim + 1))
    F_last = U[, U_index]
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


CIMLR_parameter = function (X, lambda, r, c, no.dim = NA, k = 10, cores.ratio = 1)
{
  if (is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  NITER = 30
  num = ncol(X[[1]])
  r = r
  beta = 1
  cat("Computing the multiple Kernels.\n")
  for (data_types in 1:length(X)) {
    curr_X = X[[data_types]]
    if (data_types == 1) {
      D_Kernels = multiple.kernel.cimlr(t(curr_X), cores.ratio)
    }
    else {
      D_Kernels = c(D_Kernels, multiple.kernel.cimlr(t(curr_X),
                                                     cores.ratio))
    }
  }
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
  # if (r <= 0) {
  #   r = mean(rr)
  # }
  # lambda = max(mean(rr), 0)

  lambda = lambda
  A[is.nan(A)] = 0
  S0 = max(max(distX)) - distX
  cat("Performing no network diffusion.\n")
  # S0 = network.diffusion(S0, k)
  S0 = as.matrix(S0) # new
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
    # S = network.diffusion(S, k)
    S = as.matrix(S) # new
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
        # lambda = 1.5 * lambda
        # r = r/1.01
        lambda =  lambda
        r = r
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
    # F_last = tsne(S, k = no.dim, initial_config = U[, U_index])
    F_last = U[,U_index]
  }else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U), (ncol(U) - no.dim[i] + 1))
      # F_last[i] = list(tsne(S, k = no.dim[i], initial_config = U[,U_index]))
      F_last[i] = list(U[,U_index])
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
