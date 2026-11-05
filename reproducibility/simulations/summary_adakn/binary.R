if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)

args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

i <- args[1]

# source(repro_path("annokn_source", "KF_anno.R"))
# source(repro_path("annogk_source", "GK_anno.R"))
# source(repro_path("annogk_source", "GhostKnockoff.R"))

library(AnnoKn)
library(AnnoGKR)
library(knockoff)
library(glmnet)
library(ghostbasil)

if (!requireNamespace("gam", quietly = TRUE)) {
  remotes::install_version("gam", version = "1.22-5", repos = "https://cloud.r-project.org")
}
library(gam)

source(repro_path("revision_helpers", "filter_gam.R"))
source(repro_path("revision_helpers", "filter_glm.R"))
source(repro_path("revision_helpers", "filter_EM.R"))
source(repro_path("revision_helpers", "filter_randomForest.R"))



create_gaussian_s <- function(X,Sigma,s,chol_V){
  n = nrow(X)
  p = ncol(X)
  D = diag(s)
  E = matrix(rnorm(p*n), n, p)
  P = diag(p) - solve(Sigma, D)
  V = 2*D - D%*%solve(Sigma, D)
  return(X%*%P + E%*%chol_V)
}



seed <- 1000*i

N.effect = 5000
p = 300
n0 = 30
heri = 0.05
rho = 0.5
alphalist  <- seq(0.4, 0.1, by = -0.05)
fdr_seq  <- seq(0.4, 0.1, by = -0.05)
len = length(fdr_seq)

for (heri in c(0.05, 0.1, 0.2)){
  for (p in c(300, 600, 1000)){
    set.seed(seed)
    power1 = numeric(len)
    power2 = numeric(len)
    power3 = numeric(len)
    power4 = numeric(len)
    power5 = numeric(len)
    power6 = numeric(len)
    power7 = numeric(len)
    power8 = numeric(len)
    power9 = numeric(len)
    power10 = numeric(len)
    power11 = numeric(len)

    fdr1 = numeric(len)
    fdr2 = numeric(len)
    fdr3 = numeric(len)
    fdr4 = numeric(len)
    fdr5 = numeric(len)
    fdr6 = numeric(len)
    fdr7 = numeric(len)
    fdr8 = numeric(len)
    fdr9 = numeric(len)
    fdr10 = numeric(len)
    fdr11 = numeric(len)

    amplitude = sqrt(heri/(1-heri))

    sigprob = rep(0,p)
    sigprob[1:60] = 1/(1:60)^2/(sum(1/(1:60)^2))
    # sigprob[1:60] = 1/(1:60)/(sum(1/(1:60)))
    rand = sample(1:p,n0,prob = sigprob)
    #rand = sample(1:p,n0)
    rand_sign = sample(c(-1,1), size=p, replace=TRUE)

    # AR(1) covariates
    Covariance = toeplitz(rho^(0:(p-1)))
    Sigma = Covariance
    Sigma_inv = solve(Sigma)
    s = create.solve_sdp(Sigma)
    Sigma_0 = matrix(0,2*p,2*p)
    Sigma_0[1:p,1:p] = Sigma
    Sigma_0[1:p,(p+1):(2*p)] = Sigma - diag(s)
    Sigma_0[(p+1):(2*p),1:p] = Sigma - diag(s)
    Sigma_0[(p+1):(2*p),(p+1):(2*p)] = Sigma
    Sigma_0_inv = solve(Sigma_0)
    chol_Sigma_0 = chol(Sigma_0)
    D = diag(s)
    P = diag(p) - solve(Sigma,D)
    V = 2*D - D%*%solve(Sigma,D)
    chol_V = chol(V)

    X = matrix(rnorm(N.effect*p),N.effect,p)%*%chol(Covariance)
    X = scale(X)
    beta = rep(0,p)
    beta[rand] = amplitude
    beta = beta*rand_sign
    y = X%*%beta + sqrt(n0)*rnorm(N.effect) # generate y via linear model
    print(var(X%*%beta)/var(y))
    h2e = var(X%*%beta)/var(y)
    y = (y - mean(y))/sd(y)

    X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
    X_comb = cbind(X,X_tilde)


    ############### 1. knockoff ############

    start <- Sys.time()
    set.seed(seed)
    X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
    end_check1 <- Sys.time()
    X_comb = cbind(X,X_tilde)

    mdl = cv.glmnet(X_comb,y,alpha=1)
    cvlambda = mdl$lambda.min
    beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W1 = T0-T_tilde

    end <- Sys.time()
    time1 <- as.numeric(difftime(end, start, units = "secs"))
    time1_check1 <- as.numeric(difftime(end_check1, start, units = "secs"))

    for (j in 1:len) {
      alpha = alphalist[j]
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      power1[j] = power_cal(rej1, rand)
      fdr1[j] = fdr_cal(rej1, rand)

    }


    ############### 2. AnnoKn-simple
    z <- numeric(p)
    z[1:60] = 1
    R <- scale(as.matrix(z))

    start <- Sys.time()
    set.seed(seed)
    result = AnnoKn_lite(X = X, Xk = X_tilde, y = y, R = R)
    beta = result$beta
    end <- Sys.time()
    time2 <- as.numeric(difftime(end, start, units = "secs")) + time1_check1

    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W2 = T0-T_tilde

    for (j in 1:len) {
      alpha = alphalist[j]
      tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
      rej2 = as.numeric(which(W2>=tau))
      power2[j] = power_cal(rej2, rand)
      fdr2[j] = fdr_cal(rej2, rand)
    }

    ############### 3. AnnoKn
    start <- Sys.time()
    set.seed(seed)
    result = AnnoKn(X = X, Xk = X_tilde, y = y, R = R)
    end <- Sys.time()
    beta = result$beta
    time3 <- as.numeric(difftime(end, start, units = "secs")) + time1_check1


    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W3 = T0-T_tilde

    for (j in 1:len) {

      alpha = alphalist[j]
      tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
      rej3 = as.numeric(which(W3>=tau))
      power3[j] = power_cal(rej3, rand)
      fdr3[j] = fdr_cal(rej3, rand)

    }


    ################## 4. GhostKnockoff

    ss_b = ss_se = ss_z = numeric(p)

    for(l in 1:p){
      linearmodel <- lm(y ~ X[,l])
      ss_b[l] <- coef(summary(linearmodel))[2,1]
      ss_se[l] <- coef(summary(linearmodel))[2,2]
      ss_z[l] <- coef(summary(linearmodel))[2,3]
    }

    Z       <- ss_z
    LD = cor(X)
    M = 1


    start <- Sys.time()
    set.seed(seed)
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp"
    )
    GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])

    T_0_original = GK1_lasso$T_0[[1]]
    T_k_original = GK1_lasso$T_k[[1]]

    end <- Sys.time()
    time4 <- as.numeric(difftime(end, start, units = "secs"))

    for (j in 1:len){
      threshold = alphalist[j]
      rej4 <- which(GK.filter$q <= threshold)

      power4[j] = power_cal(rej4, rand)
      fdr4[j] = fdr_cal(rej4, rand)
    }

    ############### 5. AnnoGK-simple
    start <- Sys.time()
    set.seed(seed)
    GK1_M1_anno_ps  = GK_simple(Z = Z,
                                R = R,
                                M = 1,
                                LD = LD,
                                n = N.effect,
                                ts = 'lasso')
    beta <- GK1_M1_anno_ps$beta
    T_0<-abs(beta[1:p])
    T_k<-abs(matrix(beta[-(1:p)],p,M))

    GK.filter<-GhostKnockoff.filter(T_0,T_k)

    end <- Sys.time()
    time5 <- as.numeric(difftime(end, start, units = "secs"))

    for (j in 1:len){
      threshold = alphalist[j]
      rej5 <- which(GK.filter$q <= threshold)
      power5[j] = power_cal(rej5, rand)
      fdr5[j] = fdr_cal(rej5, rand)
    }


    ############# 6. AnnoGK
    start <- Sys.time()
    set.seed(seed)
    GK1ps_anno = GK_anno_M(Z, R, M, LD, N.effect)


    beta <- GK1ps_anno$beta_final
    T_0<-GK1ps_anno$T_0
    T_k<-GK1ps_anno$T_k

    GK.filter<-GhostKnockoff.filter(T_0,T_k)

    end <- Sys.time()
    time6 <- as.numeric(difftime(end, start, units = "secs"))

    for (j in 1:len){
      threshold = alphalist[j]
      rej6 <- which(GK.filter$q <= threshold)
      power6[j] = power_cal(rej6, rand)
      fdr6[j] = fdr_cal(rej6, rand)
    }




    ############# 7. AdaKn:


    start <- Sys.time()
    set.seed(seed)
    adakn_full = filter_glm(W1, z, alpha = alphalist)

    end <- Sys.time()

    time7 <- as.numeric(difftime(end, start, units = "secs"))

    for (j in 1:len){
      rej7 <- adakn_full$rejs[[j]]
      power7[j] = power_cal(rej7, rand)
      fdr7[j] = fdr_cal(rej7, rand)
    }

    ############### 8. AdaKn (summary)

    start <- Sys.time()
    set.seed(seed)
    W_adakn <- T_0_original - as.numeric(T_k_original)

    adakn_result = filter_glm(W_adakn, z, alpha = alphalist)

    end <- Sys.time()

    time8 <- as.numeric(difftime(end, start, units = "secs"))

    for (j in 1:len){
      rej8 <- adakn_result$rejs[[j]]
      power8[j] = power_cal(rej8, rand)
      fdr8[j] = fdr_cal(rej8, rand)
    }


    result = list(power1 = power1,
                  power2 = power2,
                  power3 = power3,
                  power4 = power4,
                  power5 = power5,
                  power6 = power6,
                  power7 = power7,
                  power8 = power8,
                  fdr1 = fdr1,
                  fdr2 = fdr2,
                  fdr3 = fdr3,
                  fdr4 = fdr4,
                  fdr5 = fdr5,
                  fdr6 = fdr6,
                  fdr7 = fdr7,
                  fdr8 = fdr8,
                  time1 = time1,
                  time2 = time2,
                  time3 = time3,
                  time4 = time4,
                  time5 = time5,
                  time6 = time6,
                  time7 = time7,
                  time8 = time8,
                  h2e = h2e
    )

    path = paste0(repro_path("summary_adakn_results", "heri_final_10_binary_"),heri,"_n_",N.effect,"_p_",p, "_", i,".RData")
    # path = paste0(repro_path("annogk_results", "heri_final_10_weak_"),heri,"_n_",N.effect,"_p_",p, "_", i,".RData")

    dir <- dirname(path)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    save(result, file = path)
    print(path)

  }
}
