# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)


######################### Simulation for n > p & 1 - dim ###########################
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

i <- args[1]

####################################
# library(glmnet)
# library("xtune")
# library("SNPknock")
# suppressMessages(library("adaptMT"))
# suppressMessages(library("splines"))
# suppressMessages(library("knockoff"))
# suppressMessages(library("SNPknock"))
# suppressMessages(library("dplyr"))
# suppressMessages(library("corpcor"))
# suppressMessages(library("glmnet"))
# suppressMessages(library("MASS"))
# suppressMessages(library("gam"))
# suppressMessages(library("randomForest"))
# suppressMessages(library("tidyverse"))
# suppressMessages(library("mgcv"))
# # suppressMessages(library("flare"))
# suppressMessages(library("expm"))
# suppressMessages(library("rdist"))
# library(hdi)


source_gitfile <- function(filename){
  source(sprintf("%s.R", filename))
}
file_vec <- c(repro_path("comparison_utils", "all_other_methods"),
              repro_path("comparison_utils", "adaptive_knockoff"),
              repro_path("comparison_utils", "filter_EM"),
              repro_path("comparison_utils", "filter_gam"),
              repro_path("comparison_utils", "filter_glm"),
              repro_path("comparison_utils", "filter_randomForest"),
              repro_path("comparison_utils", "All_q_est_functions"),
              repro_path("comparison_utils", "accumulation_test_functions"))
getfile <- sapply(file_vec,source_gitfile)



source(repro_path("annokn_source", "KF_anno.R"))
#source(repro_path("archived_knockoff_source", "knockoff_anno.R"))
# source(repro_path("archived_em_source", "EM.R"))
# suppressMessages(library("knockoff"))
# suppressMessages(library("SNPknock"))
# suppressMessages(library("glmnet"))





print(paste0("######################################### iteration ", i, "  ######################################### "))

seed = 1000 * i
set.seed(seed)
n = 1000
p = 900
k = 150
d = 2

alphalist <- seq(0.3,0.05,-0.05)
Sigma = diag(rep(1,p))
amp = 3.5
binary = ''

rho = 0.5
alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)


#for(binary in c('','binary')){
for(binary in c('','binary')){
  power1 = numeric(len)
  power2 = numeric(len)
  power3 = numeric(len)
  power4 = numeric(len)
  power5 = numeric(len)
  power6 = numeric(len)
  power7 = numeric(len)
  power8 = numeric(len)
  fdr1 = numeric(len)
  fdr2 = numeric(len)
  fdr3 = numeric(len)
  fdr4 = numeric(len)
  fdr5 = numeric(len)
  fdr6 = numeric(len)
  fdr7 = numeric(len)
  fdr8 = numeric(len)
  
  # cv_matrix <- numeric(61)
  # cv_matrix_em <- numeric(61)
  
  
  Sigma = toeplitz(rho^(0:(p-1)))
  sigprob = rep(0,p)
  sigprob[1:300] = 1/(1:300)^d/(sum(1/(1:300)^d))
  nonzero = sample(1:p,k,prob = sigprob)
  beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
  y.sample = function(X) X%*%beta0+rnorm(n,0,1)
  all_res <- data.frame()
  
  ####################################
  ## HMM parameters
  ####################################
  # Number of possible states for each variable
  K=5;M=3;  
  # Marginal distribution for the first variable
  pInit = rep(1/K,K)
  # Create p-1 transition matrices
  Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
  for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
  pEmit = array(stats::runif(p*M*K),c(p,M,K))
  for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }
  
  ####################################
  ## Generating data
  ####################################
  X = sampleHMM(pInit, Q, pEmit, n=n)
  y = y.sample(X)
  
  ####################################
  ## Computing p values
  ####################################
  mdl = lm(y~X)
  pvals = summary(mdl)$coefficients[-1,4]
  
  ################ 1. AdaPT
  set.seed(seed)
  z <-1:p
  if(binary == 'binary'){
    z= ifelse(z <= 300, 0 , 1) 
  }
  pi.formulas <- paste0("ns(z, df = ", 6:10, ")") # as in Lei and Fithian
  mu.formulas <- paste0("ns(z, df = ", 6:10, ")")
  
  time1 = 0
  if(binary == ''){
    start = Sys.time()
    res <- adapt_gam(x = data.frame(z), pvals = pvals, pi_formulas = pi.formulas, mu_formulas = mu.formulas,alpha = alphalist)
    end = Sys.time()
    time1 <- as.numeric(difftime(end, start, units = "secs"))
    
    # fdr1 = c()
    # power1 = c()
    # for (j in 1:length(alphalist)){
    #   if(res$nrejs[j]>0){
    #     rej =res$rejs[[j]]
    #     fdr1 = c(fdr1,sum(beta0[rej]==0)/max(length(rej),1))
    #     power1 = c(power1,sum(beta0[rej]!=0)/max(k,1))
    #   }else{
    #     fdr1 = c(fdr1,0)
    #     power1 = c(power1,0)
    #   }
    # }
    # fdr1   <- rev(fdr1)
    # power1 <- rev(power1)
    
    for (j in 1:length(alphalist)){
      rej1 =res$rejs[[j]]
      power1[j] = power_cal(rej1, nonzero)
      fdr1[j] = fdr_cal(rej1, nonzero)
    }
    fdr1   <- rev(fdr1)
    power1 <- rev(power1)
  }
  
  ############# 2. Vanilla  knockoff
  
  start = Sys.time()
  Xk = knockoffHMM(X, pInit, Q,pEmit)
  end_check1 = Sys.time()
  mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W1 = T-T_tilde
  end_check2 = Sys.time()
  
  time2 <- as.numeric(difftime(end_check2, start, units = "secs"))
  time_check1 <- as.numeric(difftime(end_check1, start, units = "secs"))
  
  for (j in 1:len) {
    alpha = alphalist[j]
    tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
    rej2 = as.numeric(which(W1>=tau))
    power2[j] = power_cal(rej2, nonzero)
    fdr2[j] = fdr_cal(rej2, nonzero)
    
  }
  
  
  ############ 3. Adaknockoff with glm
  
  start = Sys.time()
  resj = filter_glm(W1,z,alpha =alphalist,offset=1)
  end = Sys.time()
  time3 <- as.numeric(difftime(end, start, units = "secs")) + time2
  
  for(j in 1:len){
    rej3  = resj$rejs[[j]]
    power3[j] = power_cal(rej3, nonzero)
    fdr3[j] = fdr_cal(rej3, nonzero)
  }
  
  
  ########### 4. Adaknockoff with gam
  
  start = Sys.time()
  resj = filter_gam(W1,z,alpha =alphalist,offset=1)
  end = Sys.time()
  time4 <- as.numeric(difftime(end, start, units = "secs")) + time2
  
  for(j in 1:len){
    rej4  = resj$rejs[[j]]
    power4[j] = power_cal(rej4, nonzero)
    fdr4[j] = fdr_cal(rej4, nonzero)
  }
  
  ######### 5. Adaknockoff with random forest
  start = Sys.time()
  resj = filter_randomForest(W1,z,alpha =alphalist,offset=1)
  end = Sys.time()
  time5 <- as.numeric(difftime(end, start, units = "secs")) + time2
  
  for(j in 1:len){
    rej5  = resj$rejs[[j]]
    power5[j] = power_cal(rej5, nonzero)
    fdr5[j] = fdr_cal(rej5, nonzero)
  }
  
  ######### 6. Adaknockoff with EM algorithm
  
  if(binary == ''){
    start = Sys.time()
    resj =filter_EM(W1,z,alpha =alphalist,offset=1,df = 2)
    end = Sys.time()
    time6 <- as.numeric(difftime(end, start, units = "secs")) + time2
    for(j in 1:len){
      rej6  = resj$rejs[[j]]
      power6[j] = power_cal(rej6, nonzero)
      fdr6[j] = fdr_cal(rej6, nonzero)
    }
  }else{
    start = Sys.time()
    resj =filter_EM(W1,z,alpha =alphalist,offset=1,cutoff = 0)
    end = Sys.time()
    time6 <- as.numeric(difftime(end, start, units = "secs")) + time2
    for(j in 1:len){
      rej6  = resj$rejs[[j]]
      power6[j] = power_cal(rej6, nonzero)
      fdr6[j] = fdr_cal(rej6, nonzero)
    }
  }
  
  
  ######### 7. Knockoff-anno
  set.seed(seed)
  R <- scale(as.matrix(z))
  start = Sys.time()
  result_annokn <- knockoff_anno_improved(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
  end = Sys.time()
  time7 <- as.numeric(difftime(end, start, units = "secs")) + time_check1
  beta = result_annokn$beta
  
  lambda_vals <- unlist(result_annokn$lambda_s)
  
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W3 = T-T_tilde
  
  for (j in 1:len) {
    
    alpha = alphalist[j]
    tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
    rej7 = as.numeric(which(W3>=tau))
    power7[j] = power_cal(rej7, nonzero)
    fdr7[j] = fdr_cal(rej7, nonzero)
    
  }
  
  
  ############# 8. knockoff_simple
  
  set.seed(seed)
  start <- Sys.time()
  result_annokn_lite = knockoff_simple(X = X, Xk = Xk, y = y, R = R)
  end <- Sys.time()
  time8 <- as.numeric(difftime(end, start, units = "secs")) + time_check1
  beta = result_annokn_lite$beta
  
  
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W2 = T-T_tilde
  
  for (j in 1:len) {
    
    alpha = alphalist[j]
    tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
    rej8 = as.numeric(which(W2>=tau))
    power8[j] = power_cal(rej8, nonzero)
    fdr8[j] = fdr_cal(rej8, nonzero)
    
  }
  
  
  result = list(power1 = power1, # AdaPT
                power2 = power2, # knockoff
                power3 = power3, # AdaKn (GLM)
                power4 = power4, # AdaKn (GAM) 
                power5 = power5, # AdaKn (RF) 
                power6 = power6, # AdaKn (EM)
                power7 = power7, # knockoff anno
                power8 = power8, # knockoff simple
                lambda_vals = lambda_vals,
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
                time8 = time8)
  
  
  path = paste0(repro_path("annokn_1d_results", "result_"), binary, "_amp_",amp, "_", i,".RData")
  
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  if (file.exists(path) && file.info(path)$isdir) {
    unlink(path, recursive = TRUE)
  }
  
  if (file.exists(path) && !file.info(path)$isdir) {
    file.remove(path)
  }
  
  save(result, file = path)
  
  print(path)
}