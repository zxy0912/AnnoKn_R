
suppressPackageStartupMessages({
  library(glmnet)
  library(xtune)
  library(SNPknock)
  library(adaptMT)
  library(splines)
  library(knockoff)
  library(dplyr)
  library(corpcor)
  library(MASS)
  library(gam)
  library(randomForest)
  library(tidyverse)
  library(mgcv)
  library(expm)
  library(rdist)
  library(hdi)
  library(ghostbasil)
})




obj_fun_gk_old <- function(lam, lambda_0, beta, R, n, d = 20, M = 5){
  p <- length(beta) / (M + 1)
  eta <- as.vector(R %*% lam) / d
  max_eta <- max(eta)
  w <- exp(eta - max_eta) * p / sum(exp(eta - max_eta))
  w <- w * lambda_0
  beta_mat      <- matrix(beta, nrow = p, ncol = M + 1)
  combined_beta <- rowSums(abs(beta_mat))
  part2 <- sum(w * combined_beta)
  part3 <- 0.5 * sum(lam^2) / n
  return(part2 + part3)
}


obj_fun_gk <- function(lam, lambda_0, beta, R, n, d = 20, M = 5){
  p <- length(beta) / (M + 1)
  eta <- as.vector(R %*% lam) / d
  max_eta <- max(eta)
  w <- exp(eta - max_eta) * p / sum(exp(eta - max_eta))
  beta_mat      <- matrix(beta, nrow = p, ncol = M + 1)
  combined_beta <- rowSums(abs(beta_mat))

  part1 = - (M+1) * sum(log(w)) / n
  part2 <- sum(w * lambda_0 * combined_beta)
  part3 <- 0.5 * sum(lam^2) / n
  return(part1 + part2 + part3)
}

obj_fun_gk = obj_fun_gk_old


# weight_standardized <- function(lam, R, d = 20){
#   p = nrow(R)
#   eta <- as.vector(R %*% lam) / d
#   max_eta = max(eta)
#   w <- exp(eta - max_eta) * p / (sum(exp(eta - max_eta)))
#   return(w)
# }


weight_standardized <- function(lam, R){
  print("use non-standardized weights")
  d = 10 * max(abs(R))
  p = nrow(R)
  eta <- as.vector(R %*% lam) / d
  max_eta = max(eta)
  w <- exp(eta)
  return(w)

}




# power_cal <- function(index_result, index_true){
#   if(length(index_result) == 0){
#     return(0)
#   }
#   index_common <- intersect(index_result, index_true)
#   return(length(index_common)/max(length(index_true), 1))
# }
#
#
# fdr_cal <- function(index_result, index_true){
#   if(length(index_result) == 0){
#     return(0)
#   }
#   index_common <- intersect(index_result, index_true)
#   return(1-length(index_common)/max(length(index_result), 1))
# }


################ GK-simple

GK_simple <- function(Z, R, M = 1, LD, n, ts = 'lasso'){

  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp"
  )

  GK1_lasso <- GhostKnockoff.fit(Zscore_0 = Z,
                                 N.effect = n,
                                 fit.prelim = fit.prelim,
                                 method = ts)

  result <- sol_beta_and_lam_simple_ps(fit.prelim = fit.prelim,
                                       Zscore_0 = Z,
                                       M = M,
                                       R = R,
                                       temp.A = GK1_lasso$temp.A,
                                       r_all = GK1_lasso$r_all,
                                       n = n,
                                       lambda.seq = GK1_lasso$lambda.seq,
                                       lambda_0 = GK1_lasso$lambda)

  return(result)

}


sol_beta_and_lam_simple_ps = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, n,
                                      lambda.seq, lambda_0, init_scale = NULL, maxiter = 100,
                                      verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0
  if(is.null(init_scale)){
    init_scale = 0
  }

  p = length(Zscore_0)
  Prior = numeric((M+1)*p)+1

  fit_basil <- ghostbasil(
    temp.A, r_all,
    user.lambdas        = lambda.seq* init_scale,
    alpha               = 1,
    penalty= Prior,
    delta.strong.size   = max(1, min(500, length(r_all)/20)),
    max.strong.size     = length(r_all),
    n.threads           = 1,
    use.strong.rule     = FALSE
  )

  beta <- fit_basil$betas[, ncol(fit_basil$betas)]


  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100


  lambda_s_all <- numeric()
  objs = c()

  init_lam <- numeric(ncol(R)) + 5

  for(t in 1:100){

    res <- optim(
      par   = init_lam,
      fn    = obj_fun_gk,
      method = "BFGS",
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n,
      M = M
    )

    lambda_s = res$par
    print(lambda_s)
    lambda_s_all <- append(lambda_s_all, lambda_s)

    w = weight_standardized(lambda_s, R)

    Prior = rep(w,M+1)

    fit_basil <- ghostbasil(
      temp.A, r_all,
      user.lambdas        = lambda.seq* sum(Prior) / ((M+1)*p),
      alpha               = 1,
      penalty= Prior,
      delta.strong.size   = max(1, min(500, length(r_all)/20)),
      max.strong.size     = length(r_all),
      n.threads           = 1,
      use.strong.rule     = FALSE
    )

    beta      <- fit_basil$betas[, ncol(fit_basil$betas)]


    obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
    obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)


    cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")

    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }

    lambda_s_old = lambda_s

    beta_old = beta
  }
  n.G<-length(Zscore_0)
  P.each<-fit.prelim$P.each
  Normal_50Studies<-fit.prelim$Normal_50Studies
  A<-as.matrix(fit.prelim$A)
  A.left<-fit.prelim$A.left
  V.left<-fit.prelim$V.left

  T_0<-list();T_k<-list()

  Normal_k<-matrix(V.left%*%matrix(rnorm(ncol(V.left)),ncol(V.left),1),nrow=n.G)
  GK.Zscore_0<-Zscore_0
  GK.Zscore_k<-as.vector(P.each%*%GK.Zscore_0)+Normal_k

  r<-GK.Zscore_0/sqrt(n)#sqrt(n-1+GK.Zscore_0^2)
  r_k<-as.vector(GK.Zscore_k/sqrt(n))#sqrt(n-1+GK.Zscore_k^2))
  r_all<-as.matrix(c(r,r_k))

  nfold<-5
  nA<-n*(nfold-1)/nfold;nB<-n/nfold
  temp.left<-sqrt(nB/nA/n)*as.matrix(A.left)
  r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  r_all_B<-(r_all*n-r_all_A*nA)/nB
  shrink=0.01#seq(0.05,1,0.05)
  beta.all<-c();parameter.set<-c()
  k<-1
  #temp.A<-(1-shrink[k])*A+diag(shrink[k],nrow(A))
  temp.A<-A+diag(shrink[k],nrow(A))
  fit.basil<-try(ghostbasil(temp.A, r_all_A, alpha=1, penalty= Prior, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=FALSE),silent=T)
  parameter.set<-rbind(parameter.set,cbind(fit.basil$lmdas,shrink[k]))
  beta.all<-cbind(beta.all,fit.basil$betas)

  Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
  f.lambda<-apply(beta.all,2,Get.f)
  f.lambda[is.na(f.lambda)]<--Inf
  #beta<-beta.all[,which.max(f.lambda)]
  parameter<-parameter.set[which.max(f.lambda),]
  temp.A<-(1-parameter[2])*A+diag(parameter[2],nrow(A))

  lambda.all<-fit.basil$lmdas
  lambda<-fit.basil$lmdas[which.max(f.lambda)]
  lambda.seq <- lambda.all[lambda.all > lambda]
  lambda.seq <- c(lambda.seq, lambda)

  fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, penalty= Prior, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=FALSE)
  beta<-fit.basil$betas[,ncol(fit.basil$betas)]

  list(beta = beta, objs = objs, lambda_s = lambda_s, Prior = Prior)
}









############################## AnnoGK ##############################


sol_beta_and_lam_ps = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, n,
                               lambda.seq, lambda_0, init_scale = NULL, maxiter = 100,
                               verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0
  if(is.null(init_scale)){
    init_scale = 0
  }

  p = length(Zscore_0)
  Prior = numeric((M+1)*p)+1

  fit_basil <- ghostbasil(
    temp.A, r_all,
    user.lambdas        = lambda.seq* init_scale,
    alpha               = 1,
    penalty= Prior,
    delta.strong.size   = max(1, min(500, length(r_all)/20)),
    max.strong.size     = length(r_all),
    n.threads           = 1,
    use.strong.rule     = FALSE
  )

  beta <- fit_basil$betas[, ncol(fit_basil$betas)]


  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100


  lambda_s_all <- numeric()
  objs = c()

  init_lam <- numeric(ncol(R)) + 5

  for(t in 1:100){

    res <- optim(
      par   = init_lam,
      fn    = obj_fun_gk,
      method = "BFGS",
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n,
      M = M
    )

    lambda_s = res$par
    print(lambda_s)
    lambda_s_all <- append(lambda_s_all, lambda_s)

    w = weight_standardized(lambda_s, R)

    Prior = rep(w,M+1)

    fit_basil <- ghostbasil(
      temp.A, r_all,
      user.lambdas        = lambda.seq* sum(Prior) / ((M+1)*p),
      alpha               = 1,
      penalty= Prior,
      delta.strong.size   = max(1, min(500, length(r_all)/20)),
      max.strong.size     = length(r_all),
      n.threads           = 1,
      use.strong.rule     = FALSE
    )

    beta      <- fit_basil$betas[, ncol(fit_basil$betas)]


    obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
    obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)


    cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")

    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }

    lambda_s_old = lambda_s

    beta_old = beta
  }

  list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior, lambda_s = lambda_s)
}


############### objective function used for AnnoGK ##############################

cv_func <- function(r_all, temp.A, A.left, lambda.seq, Prior, n, nfold = 5){
  nA<-n*(nfold-1)/nfold;nB<-n/nfold
  temp.left<-sqrt(nB/nA/n)*as.matrix(A.left)
  r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  r_all_B<-(r_all*n-r_all_A*nA)/nB

  fit.basil<-ghostbasil(temp.A, r_all_A, user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)),
                        penalty= Prior, max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
  beta<-fit.basil$betas[,ncol(fit.basil$betas)]
  # print(hist(beta))
  if(all(beta == 0) == TRUE){
    return(-Inf)
  }
  Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
  return(Get.f(beta))
}




# GK_anno <- function(Z, R, M = 1, LD, n, ts = 'lasso'){
#
#   p = nrow(LD)
#   var_para <- apply(R, 2, var)
#
#   if(all(is.na(R)) | sum(var_para!=0) == 0){
#     print("no informative annotations")
#     fit.prelim <- GhostKnockoff.prelim(
#       cor.G   = LD,
#       M       = M,
#       method  = "sdp"
#     )
#     GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method='lasso')
#     # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
#     T_0 = GK1_lasso$T_0[[1]]
#     T_k = GK1_lasso$T_k[[1]]
#
#     return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
#                 obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
#   }
#
#   R <- R[, var_para!=0]
#   R = scale(R)
#
#   # set.seed(123)
#
#   fit.prelim <- GhostKnockoff.prelim(
#     cor.G   = LD,
#     M       = M,
#     method  = "sdp"
#   )
#
#   GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method=ts)
#   temp.A = GK1_lasso$temp.A
#   r_all = GK1_lasso$r_all
#   lambda.seq = GK1_lasso$lambda.seq
#   lambda = GK1_lasso$lambda
#
#   coeff_list <- exp(seq(-2, 2, by = 0.2))
#
#   Prior_final = numeric(2 * p) + 1
#   CV_BEST = -Inf
#   lambda_BEST = 0
#   lambda_s_final = NA
#   beta_final = NA
#   obj_cor_list <- numeric()
#
#   for(coeff in coeff_list){
#
#     GKanno_M1_anno_ps <- sol_beta_and_lam_ps(fit.prelim = fit.prelim,
#                                              Zscore_0 = Z,
#                                              M = M,
#                                              R = R,
#                                              temp.A = temp.A,
#                                              r_all = r_all,
#                                              n = n,
#                                              lambda.seq = coeff*lambda.seq,
#                                              lambda_0 = coeff*lambda)
#
#     obj_cor <- cv_func(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq,
#                        Prior = GKanno_M1_anno_ps$Prior, n = n, nfold = 5)
#
#     obj_cor_list <- append(obj_cor_list, obj_cor)
#
#     if(obj_cor > CV_BEST){
#       lambda_BEST = coeff*lambda
#       CV_BEST = obj_cor
#       Prior_final = GKanno_M1_anno_ps$Prior
#       lambda_s_final = GKanno_M1_anno_ps$lambda_s
#       beta_final = GKanno_M1_anno_ps$beta
#     }
#   }
#
#   T_0<-abs(beta_final[1:p])
#   T_k<-abs(matrix(beta_final[-(1:p)],p,M))
#
#   return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
#               obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
# }



################# main function for AnnoGK, allowing multiple knockoff copies ##########

GK_anno_M <- function(Z, R, M = 1, LD, n, ts = 'lasso'){

  p = nrow(LD)
  var_para <- apply(R, 2, var)

  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp"
    )
    GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method='lasso')
    # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    T_0 = GK1_lasso$T_0[[1]]
    T_k = GK1_lasso$T_k[[1]]

    return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
                obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
  }

  R <- R[, var_para!=0]
  R = scale(R)

  # set.seed(123)

  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp"
  )

  GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method=ts)
  temp.A = GK1_lasso$temp.A
  r_all = GK1_lasso$r_all
  lambda.seq = GK1_lasso$lambda.seq
  lambda = GK1_lasso$lambda

  coeff_list <- exp(seq(-2, 2, by = 0.2))

  Prior_final = numeric((M+1) * p) + 1
  CV_BEST = -Inf
  lambda_BEST = 0
  lambda_s_final = NA
  beta_final = NA
  obj_cor_list <- numeric()

  for(coeff in coeff_list){

    GKanno_M1_anno_ps <- sol_beta_and_lam_ps(fit.prelim = fit.prelim,
                                             Zscore_0 = Z,
                                             M = M,
                                             R = R,
                                             temp.A = temp.A,
                                             r_all = r_all,
                                             n = n,
                                             lambda.seq = coeff*lambda.seq,
                                             lambda_0 = coeff*lambda)

    obj_cor <- cv_func(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq,
                       Prior = GKanno_M1_anno_ps$Prior, n = n, nfold = 5)

    obj_cor_list <- append(obj_cor_list, obj_cor)

    if(obj_cor > CV_BEST){
      lambda_BEST = coeff*lambda
      CV_BEST = obj_cor
      Prior_final = GKanno_M1_anno_ps$Prior
      lambda_s_final = GKanno_M1_anno_ps$lambda_s
      beta_final = GKanno_M1_anno_ps$beta
    }
  }

  T_0<-abs(beta_final[1:p])
  T_k<-abs(matrix(beta_final[-(1:p)],p,M))

  return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
              obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
}







################# modify GK_anno_M to use different effect sizes to modify r,
################# but omit different sample sizes to modify correlation matrix

sol_beta_and_lam_ps_N = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, n,
                               lambda.seq, lambda_0, init_scale = NULL, maxiter = 100,
                               verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0
  if(is.null(init_scale)){
    init_scale = 0
  }

  p = length(Zscore_0)
  Prior = numeric((M+1)*p)+1

  fit_basil <- ghostbasil(
    temp.A, r_all,
    user.lambdas        = lambda.seq* init_scale,
    alpha               = 1,
    penalty= Prior,
    delta.strong.size   = max(1, min(500, length(r_all)/20)),
    max.strong.size     = length(r_all),
    n.threads           = 1,
    use.strong.rule     = FALSE
  )

  beta <- fit_basil$betas[, ncol(fit_basil$betas)]


  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100


  lambda_s_all <- numeric()
  objs = c()

  init_lam <- numeric(ncol(R)) + 5

  for(t in 1:100){

    res <- optim(
      par   = init_lam,
      fn    = obj_fun_gk,
      method = "BFGS",
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n,
      M = M
    )

    lambda_s = res$par
    print(lambda_s)
    lambda_s_all <- append(lambda_s_all, lambda_s)

    w = weight_standardized(lambda_s, R)

    Prior = rep(w,M+1)

    fit_basil <- ghostbasil(
      temp.A, r_all,
      user.lambdas        = lambda.seq* sum(Prior) / ((M+1)*p),
      alpha               = 1,
      penalty= Prior,
      delta.strong.size   = max(1, min(500, length(r_all)/20)),
      max.strong.size     = length(r_all),
      n.threads           = 1,
      use.strong.rule     = FALSE
    )

    beta      <- fit_basil$betas[, ncol(fit_basil$betas)]


    obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
    obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)


    cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")

    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }

    lambda_s_old = lambda_s

    beta_old = beta
  }

  list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior, lambda_s = lambda_s)
}


cv_func_N <- function(r_all, temp.A, A.left, lambda.seq, Prior, N.effect, nfold = 5){
  nA<-N.effect*(nfold-1)/nfold;nB<-N.effect/nfold
  temp.left<-sqrt(nB/nA/N.effect)*as.matrix(A.left)
  r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  r_all_B<-(r_all*N.effect-r_all_A*nA)/nB

  fit.basil<-ghostbasil(temp.A, r_all_A, user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)),
                        penalty= Prior, max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
  beta<-fit.basil$betas[,ncol(fit.basil$betas)]
  # print(hist(beta))
  if(all(beta == 0) == TRUE){
    return(-Inf)
  }
  Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
  return(Get.f(beta))
}


GK_anno_M_N <- function(Z, R, M = 1, LD, N.effect, ts = 'lasso'){

  n = median(N.effect)
  p = nrow(LD)
  var_para <- apply(R, 2, var)

  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp"
    )
    GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    T_0 = GK1_lasso$T_0[[1]]
    T_k = GK1_lasso$T_k[[1]]

    return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
                obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
  }

  R <- R[, var_para!=0]
  R = scale(R)

  # set.seed(123)

  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp"
  )

  GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method=ts)
  temp.A = GK1_lasso$temp.A
  r_all = GK1_lasso$r_all
  lambda.seq = GK1_lasso$lambda.seq
  lambda = GK1_lasso$lambda

  coeff_list <- exp(seq(-2, 2, by = 0.2))

  Prior_final = numeric((M+1) * p) + 1
  CV_BEST = -Inf
  lambda_BEST = 0
  lambda_s_final = NA
  beta_final = NA
  obj_cor_list <- numeric()

  for(coeff in coeff_list){

    GKanno_M1_anno_ps <- sol_beta_and_lam_ps_N(fit.prelim = fit.prelim,
                                             Zscore_0 = Z,
                                             M = M,
                                             R = R,
                                             temp.A = temp.A,
                                             r_all = r_all,
                                             n = n,
                                             lambda.seq = coeff*lambda.seq,
                                             lambda_0 = coeff*lambda)

    obj_cor <- cv_func_N(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq,
                       Prior = GKanno_M1_anno_ps$Prior, N.effect = N.effect, nfold = 5)

    obj_cor_list <- append(obj_cor_list, obj_cor)

    if(obj_cor > CV_BEST){
      lambda_BEST = coeff*lambda
      CV_BEST = obj_cor
      Prior_final = GKanno_M1_anno_ps$Prior
      lambda_s_final = GKanno_M1_anno_ps$lambda_s
      beta_final = GKanno_M1_anno_ps$beta
    }
  }

  T_0<-abs(beta_final[1:p])
  T_k<-abs(matrix(beta_final[-(1:p)],p,M))

  return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
              obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
}







################# considering sample size for different SNPs





sol_beta_and_lam_ps_dss = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, N.effect,
                               lambda.seq, lambda_0, init_scale = NULL, maxiter = 100,
                               verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0
  if(is.null(init_scale)){
    init_scale = 0
  }
  n = median(N.effect)

  p = length(Zscore_0)
  Prior = numeric((M+1)*p)+1

  fit_basil <- ghostbasil(
    temp.A, r_all,
    user.lambdas        = lambda.seq* init_scale,
    alpha               = 1,
    penalty= Prior,
    delta.strong.size   = max(1, min(500, length(r_all)/20)),
    max.strong.size     = length(r_all),
    n.threads           = 1,
    use.strong.rule     = FALSE
  )

  beta <- fit_basil$betas[, ncol(fit_basil$betas)]


  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100


  lambda_s_all <- numeric()
  objs = c()

  init_lam <- numeric(ncol(R)) + 5

  for(t in 1:100){

    res <- optim(
      par   = init_lam,
      fn    = obj_fun_gk,
      method = "BFGS",
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n,
      M = M
    )

    lambda_s = res$par
    print(lambda_s)
    lambda_s_all <- append(lambda_s_all, lambda_s)

    w = weight_standardized(lambda_s, R)

    Prior = rep(w,M+1)

    fit_basil <- ghostbasil(
      temp.A, r_all,
      user.lambdas        = lambda.seq* sum(Prior) / ((M+1)*p),
      alpha               = 1,
      penalty= Prior,
      delta.strong.size   = max(1, min(500, length(r_all)/20)),
      max.strong.size     = length(r_all),
      n.threads           = 1,
      use.strong.rule     = FALSE
    )

    beta      <- fit_basil$betas[, ncol(fit_basil$betas)]


    obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
    obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)


    cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")

    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }

    lambda_s_old = lambda_s

    beta_old = beta
  }

  list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior, lambda_s = lambda_s)
}




GK_anno_dss <- function(Z, R, M = 1, LD, N.effect, ts = 'lasso'){

  p = nrow(LD)

  if(length(N.effect) == 1){
    N.effect = rep(N.effect, p)
  }

  var_para <- apply(R, 2, var)
  N.median = median(N.effect)

  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp"
    )
    GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    T_0 = GK1_lasso$T_0[[1]]
    T_k = GK1_lasso$T_k[[1]]

    return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
                obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
  }

  R <- R[, var_para!=0]
  R = scale(R)

  # set.seed(123)

  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp"
  )

  GK1_lasso <- GhostKnockoff.lambda(Z, N.effect, fit.prelim, method=ts)
  A = GK1_lasso$A
  modi_para <- GK1_lasso$modi_para
  r_all = GK1_lasso$r_all
  lambda.seq = GK1_lasso$lambda.seq
  lambda = GK1_lasso$lambda


  #### check point
  # GK1_lasso <- GhostKnockoff.fit(Z, mean(N.effect), fit.prelim, method=ts)
  # temp.A = GK1_lasso$temp.A
  # r_all = GK1_lasso$r_all
  # lambda.seq = GK1_lasso$lambda.seq
  # lambda = GK1_lasso$lambda


  ################# transform temp.A and r_all:
  N.matrix = diag(rep(N.effect, (M+1))/N.median)
  N.matrix.sqrt <- diag(sqrt(rep(N.effect, M+1)/N.median))
  A.modi = N.matrix.sqrt %*% A %*% N.matrix.sqrt
  r_all.modi = N.matrix %*% r_all

  temp.A.modi<-(1-modi_para)*A.modi+diag(modi_para,nrow(A.modi))
  ###################################################


  coeff_list <- exp(seq(-2, 2, by = 0.2))

  Prior_final = numeric((M+1) * p) + 1
  CV_BEST = -Inf
  lambda_BEST = 0
  lambda_s_final = NA
  beta_final = NA
  obj_cor_list <- numeric()

  for(coeff in coeff_list){

    # GKanno_M1_anno_ps <- sol_beta_and_lam_ps(fit.prelim = fit.prelim,
    #                                          Zscore_0 = Z,
    #                                          M = M,
    #                                          R = R,
    #                                          temp.A = temp.A,
    #                                          r_all = r_all,
    #                                          n = n,
    #                                          lambda.seq = coeff*lambda.seq,
    #                                          lambda_0 = coeff*lambda)

    GKanno_M1_anno_ps <- sol_beta_and_lam_ps_dss(fit.prelim = fit.prelim,
                                                 Zscore_0 = Z,
                                                 M = M,
                                                 R = R,
                                                 temp.A = temp.A.modi,
                                                 r_all = r_all.modi,
                                                 N.effect = N.effect,
                                                 lambda.seq = coeff*lambda.seq,
                                                 lambda_0 = coeff*lambda)


    # obj_cor_old <- cv_func(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq,
    #                    Prior = GKanno_M1_anno_ps$Prior, n = n, nfold = 5)

    obj_cor <- cv_func_dss(r_all = r_all, temp.A.modi = temp.A.modi, A.left = fit.prelim$A.left,
                           lambda.seq = coeff*lambda.seq,
                           Prior = GKanno_M1_anno_ps$Prior, N.effect = N.effect,
                           M = M, nfold = 5)

    obj_cor_list <- append(obj_cor_list, obj_cor)

    if(obj_cor > CV_BEST){
      lambda_BEST = coeff*lambda
      CV_BEST = obj_cor
      Prior_final = GKanno_M1_anno_ps$Prior
      lambda_s_final = GKanno_M1_anno_ps$lambda_s
      beta_final = GKanno_M1_anno_ps$beta
    }
  }

  T_0<-abs(beta_final[1:p])
  T_k<-abs(matrix(beta_final[-(1:p)],p,M))

  return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
              obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
}






cv_func_dss <- function(r_all, temp.A.modi, A.left, lambda.seq, Prior, N.effect, M = M, nfold = 5){
  # nA<-n*(nfold-1)/nfold;nB<-n/nfold
  # temp.left<-sqrt(nB/nA/n)*as.matrix(A.left)
  # r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  # r_all_B<-(r_all*n-r_all_A*nA)/nB

  N.median = median(N.effect)
  nA<-N.effect*(nfold-1)/nfold;nB<-N.effect/nfold
  temp.left<-sqrt(nB/nA/N.effect)*as.matrix(A.left)
  r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  r_all_B<-(r_all*N.effect-r_all_A*nA)/nB

  ###################################
  N.matrix = diag(rep(N.effect, (M+1))/N.median)
  # N.matrix.sqrt <- diag(sqrt(rep(N.effect, M+1)/N.median))
  # temp.A.modi = N.matrix.sqrt %*% temp.A %*% N.matrix.sqrt
  r_all_A.modi = N.matrix %*% r_all_A
  ####################################

  fit.basil<-ghostbasil(temp.A.modi, r_all_A.modi, user.lambdas=lambda.seq, alpha=1,
                        delta.strong.size = max(1,min(500,length(r_all_A)/20)),
                        penalty= Prior, max.strong.size = nrow(temp.A.modi),n.threads=1,use.strong.rule=F)
  beta<-fit.basil$betas[,ncol(fit.basil$betas)]
  # print(hist(beta))
  if(all(beta == 0) == TRUE){
    return(-Inf)
  }
  Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A.modi%*%x))}
  return(Get.f(beta))
}










GhostKnockoff.lambda<-function(Zscore_0, N.effect, fit.prelim, method='lasso',type='fdr',M.fwer=50){
  Zscore_0<-as.matrix(Zscore_0)
  N.effect<-as.vector(N.effect)
  N.median <- median(N.effect)
  Zscore_0[is.na(Zscore_0)]<-0
  N.effect[is.na(N.effect)]<-Inf

  M<-fit.prelim$M
  n.G<-nrow(Zscore_0)
  P.each<-fit.prelim$P.each
  Normal_50Studies<-fit.prelim$Normal_50Studies
  A<-as.matrix(fit.prelim$A)
  A.left<-fit.prelim$A.left
  V.left<-fit.prelim$V.left

  if(type=='fdr'){M.rep<-1}
  if(type=='fwer'){M.rep<-M.fwer}

  # T_0<-list();T_k<-list()
  # kappa<-list();tau<-list()
  for(m in 1:M.rep){
    #Normal_k<-matrix(Normal_50Studies[,m],nrow=n.G)
    Normal_k<-matrix(V.left%*%matrix(rnorm(ncol(V.left)),ncol(V.left),1),nrow=n.G)
    GK.Zscore_0<-Zscore_0
    GK.Zscore_k<-as.vector(P.each%*%GK.Zscore_0)+Normal_k
    #
    # if(method=='marginal'){
    #   T_0[[m]]<-(GK.Zscore_0)^2
    #   T_k[[m]]<-(GK.Zscore_k)^2
    # }
    if(method=='lasso'){
      #calculate importance score
      r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
      r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
      r_all<-as.matrix(c(r,r_k))

      nfold<-5
      nA<-N.effect*(nfold-1)/nfold;nB<-N.effect/nfold
      temp.left<-sqrt(nB/nA/N.effect)*as.matrix(A.left)
      r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
      r_all_B<-(r_all*N.effect-r_all_A*nA)/nB
      shrink=0.01#seq(0.05,1,0.05)
      beta.all<-c();parameter.set<-c()
      k<-1
      #temp.A<-A+diag(shrink[k],nrow(A))

      ####################### modify the matrix, although we are considering only the training set, the 0.8
      ####################### is removed from both the numerator and denominator
      N.matrix = diag(rep(N.effect, (M+1))/N.median)
      N.matrix.sqrt <- diag(sqrt(rep(N.effect, M+1)/N.median))
      A.modi = N.matrix.sqrt %*% A %*% N.matrix.sqrt
      r_all_A.modi = N.matrix %*% r_all_A

      temp.A.modi<-A.modi+diag(shrink[k],nrow(A))
      #######################

      fit.basil<-try(ghostbasil(temp.A.modi, r_all_A.modi, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A.modi),n.threads=1,use.strong.rule=FALSE),silent=T)
      parameter.set<-rbind(parameter.set,cbind(fit.basil$lmdas,shrink[k]))
      beta.all<-cbind(beta.all,fit.basil$betas)

      Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A.modi%*%x))}
      f.lambda<-apply(beta.all,2,Get.f)
      f.lambda[is.na(f.lambda)]<--Inf
      #beta<-beta.all[,which.max(f.lambda)]
      parameter<-parameter.set[which.max(f.lambda),]
      temp.A<-(1-parameter[2])*A+diag(parameter[2],nrow(A))

      lambda.all<-fit.basil$lmdas
      lambda<-fit.basil$lmdas[which.max(f.lambda)]
      lambda.seq <- lambda.all[lambda.all > lambda]
      lambda.seq <- c(lambda.seq, lambda)

      # fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=FALSE)
      # beta<-fit.basil$betas[,ncol(fit.basil$betas)]
      #
      # T_0[[m]]<-abs(beta[1:n.G])
      # T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    # if(method=='lasso.approx.lambda'){
    #   #calculate importance score
    #   r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
    #   r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
    #   r_all<-as.matrix(c(r,r_k))
    #
    #   N_eff <- as.numeric(N.effect[1])
    #   lambda_max<-max(abs(rnorm(length(r_all))))/sqrt(as.numeric(N_eff))
    #   epsilon <- .0001
    #   K <- 100
    #   lambda.all <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
    #                               length.out = K)), digits = 10)
    #   lambda<-lambda_max*0.6
    #   lambda.seq <- lambda.all[lambda.all > lambda]
    #   lambda.seq <- c(lambda.seq, lambda)
    #
    #   temp.A<-A+0.01*diag(1,nrow(A))
    #   fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
    #   beta<-fit.basil$betas[,ncol(fit.basil$betas)]
    #
    #   T_0[[m]]<-abs(beta[1:n.G])
    #   T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    # }
    # if(method=='susie'){
    #   fitted_rss <- suppressMessages(susie_rss(z=c(GK.Zscore_0,as.vector(GK.Zscore_k)), R=as.matrix(A), L = min(10,length(GK.Zscore_0))))
    #   fitted_vars<-summary(fitted_rss)$vars
    #   beta<-fitted_vars[order(fitted_vars[,1]),2]#*c(GK.Zscore_0,as.vector(GK.Zscore_k))^2
    #   T_0[[m]]<-abs(beta[1:n.G])
    #   T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    # }
    # MK.stat<-MK.statistic(T_0[[m]],T_k[[m]])
    # kappa[[m]]<-MK.stat[,'kappa']
    # tau[[m]]<-MK.stat[,'tau']
  }

  return(list(A=A,
              modi_para = parameter[2],
              r_all=r_all,
              lambda.seq=lambda.seq,
              lambda = lambda))
}
