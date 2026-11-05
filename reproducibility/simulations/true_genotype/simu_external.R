# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)



args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

iter <- args[1]

seed = 1234*iter

source(repro_path("annokn_source", "KF_anno.R"))
source(repro_path("annogk_source", "GK_anno.R"))
source(repro_path("annogk_source", "GhostKnockoff.R"))
packageVersion("Matrix") 

regularization = ''
regularization = 'regu'


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





library(tidyr)
library(dplyr)
library(genio)

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


maf_cal <- function(x){
  freq <- sum(x, na.rm = TRUE)/(2 * sum(!is.na(x)))
  maf <- ifelse(freq < 0.5, freq, 1 - freq)
  return(maf)
}



N.effect = 10000 # sample size
n1 = 1000 # sample size used for clustering
# p = 200 # number of variable
h2e = 0.1 # heritability
alphae = 0.1 # proportion of causal SNPs
p1 = 100 # number of risk SNPs
# rho = 0.5 # maximum correlation between SNPs
M = 1 # number of knockoff copies
alphalist  <- seq(0.4, 0.05, by = -0.05)
fdr_seq  <- seq(0.4, 0.05, by = -0.05)
len = length(fdr_seq)


#### load the information for each risk region
ancestry = 'EUR'
chrid = '1'
path <- paste0(repro_path("ukbb", ""),ancestry,"_0.05maf/bychr/chr_",chrid,"_2w_common")
genotype <- read_plink(path)

snp_number = read.table(file = paste0(repro_path("ukbb", "EUR_0.05maf/real_simulation/region_information/chr"), chrid,"_common.txt"))
good_iter <- which(snp_number$V3 > 100)

ld_block <- read.table(repro_path("ukbb", "EUR/fourier_ls-all.bed"), header = TRUE)
ld_block = ld_block[ld_block$chr == paste0("chr", chrid),]



path <- paste(repro_path("reference_plink", "1KG_chr"), chrid,sep='')
ref_chr <- read_plink(path)





######### start the iteration simulations:

for(h2e in c(0.05, 0.1, 0.2)){
  
  for(alphae in c(0.1, 0.2, 0.3)){
    set.seed(seed)
    i = sample(good_iter, 1)
    
    start = ld_block[i,2]
    end = ld_block[i,3]
    index = which(genotype$bim$pos >= start & genotype$bim$pos <= end)
    
    X = t(genotype$X[index, 1:N.effect])
    bim_ukb <- genotype$bim[index,]
    maf <- apply(X, 2, maf_cal)
    X = bedNA(X)
    
    #### use hierarchical clustering to cluster SNPs
    
    R = cor(X[1:n1,])
    R2 <- R^2
    d <- as.dist(1 - R2)
    hc <- hclust(d, method = "average")
    # plot(hc, labels = FALSE, main = "SNP clustering by LD")
    clusters <- cutree(hc, h = 0.75)  # e.g., cluster SNPs with r² > 0.6
    length(table(clusters))
    mean(table(clusters)) # average size for each cluster
    
    
    
    ######### select representative SNP
    
    df <- data.frame(
      SNP = names(clusters),
      maf = maf,
      cluster = clusters,
      index = 1:length(maf)
    )
    
    representatives <- df %>%
      group_by(cluster) %>%
      slice_max(order_by = maf, n = 1, with_ties = FALSE)
    
    head(representatives)
    dim(representatives)
    index <- representatives$index[order(representatives$index)]
    X = X[, index]
    bim_ukb <- bim_ukb[index, ]
    
    
    ##########
    
    
    X = scale(X)
    p = nrow(representatives)
    betae <- rep(0,p)
    # rand <- sample(p,alphae*p)
    
    sigprob = rep(0,p)
    sigprob[1:p1] = 1/(1:p1)^2/(sum(1/(1:p1)^2))
    rand = sample(1:p,alphae*p,prob = sigprob)
    
    print(rand)
    betae[rand] <- rnorm(length(rand),0,sqrt(h2e/p/alphae))
    y <- X%*%betae+rnorm(N.effect,0,sqrt(1-h2e))
    y = (y - mean(y))/sd(y)
    
    
    ss_b = ss_se = ss_z = numeric(p)
    
    for(l in 1:p){
      linearmodel <- lm(y ~ X[,l])
      ss_b[l] <- coef(summary(linearmodel))[2,1]
      ss_se[l] <- coef(summary(linearmodel))[2,2]
      ss_z[l] <- coef(summary(linearmodel))[2,3]
    }
    
    Z       <- ss_z  
    LD <- cor(X)
    n0 = length(rand)
    
    
    ##################
    
    power1 = numeric(len) ## knockoff
    power2 = numeric(len) ## knockoff-simple
    power3 = numeric(len) ## knockoff-anno
    power4 = numeric(len) ## Ghostknockoff 
    power5 = numeric(len) ## GK-simple pseudo sum
    power6 = numeric(len) ## GK-anno pseudo sum
    power7 = numeric(len) ## AdaKn + EM
    power8 = numeric(len) ## AdaKn + random forest
    
    
    fdr1 = numeric(len) ## knockoff
    fdr2 = numeric(len) ## knockoff-simple
    fdr3 = numeric(len) ## knockoff-anno
    fdr4 = numeric(len) ## Ghostknockoff 
    fdr5 = numeric(len) ## GK-simple pseudo sum
    fdr6 = numeric(len) ## GK-anno pseudo sum
    fdr7 = numeric(len) ## AdaKn + EM
    fdr8 = numeric(len) ## AdaKn + random forest
    
    
    ############### 1. knockoff ############
    start <- Sys.time()
    set.seed(seed)
    X_tilde = create.gaussian(X, rep(0,p), LD) # generate knockoff variable
    X_comb = cbind(X,X_tilde)
    
    mdl = cv.glmnet(X_comb,y,alpha=1)
    cvlambda = mdl$lambda.min
    beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W1 = T0-T_tilde
    
    end <- Sys.time()
    
    for (j in 1:len) {
      alpha = alphalist[j]
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      power1[j] = power_cal(rej1, rand)
      fdr1[j] = fdr_cal(rej1, rand)
      
    }
    time1 <- as.numeric(difftime(end, start, units = "secs"))
    cat("knockoff time:", time1, "seconds\n")
    
    
    
    
    ############### 2. knockoff simple with annotation
    
    z <-1:p
    R <- scale(as.matrix(z))
    
    start <- Sys.time()
    set.seed(seed)
    result = knockoff_simple(X = X, Xk = X_tilde, y = y, R = R)
    beta = result$beta
    
    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W2 = T0-T_tilde
    
    end <- Sys.time()
    time2 <- as.numeric(difftime(end, start, units = "secs"))
    cat("knockoff-simple time:", time2, "seconds\n")
    
    for (j in 1:len) {
      alpha = alphalist[j]
      tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
      rej2 = as.numeric(which(W2>=tau))
      power2[j] = power_cal(rej2, rand)
      fdr2[j] = fdr_cal(rej2, rand)
    }
    
    ############### 3. knockoff-anno with annotation
    
    start <- Sys.time()
    set.seed(seed)
    result = knockoff_anno_improved(X = X, Xk = X_tilde, y = y, R = R)
    beta = result$beta
    
    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W3 = T0-T_tilde
    
    
    end <- Sys.time()
    time3 <- as.numeric(difftime(end, start, units = "secs"))
    cat("knockoff-anno time:", time3, "seconds\n")
    
    
    for (j in 1:len) {
      
      alpha = alphalist[j]
      tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
      rej3 = as.numeric(which(W3>=tau))
      power3[j] = power_cal(rej3, rand)
      fdr3[j] = fdr_cal(rej3, rand)
      
    }
    
    
    
    
    ########### external reference panel:
    
    
    
    snp_name <- colnames(X)
    snp_index <- match(snp_name, ref_chr$bim$id)
    X_ref <- ref_chr$X[snp_index,]
    bim_ref <- ref_chr$bim[snp_index,]
    X_ref <- bedNA(t(X_ref))
    X_ref <- scale(X_ref)
    
    LD_ref <- cor(X_ref)
    
    
    library(susieR)
    lambda = estimate_s_rss(Z, LD_ref, n = N.effect)
    print(lambda)
    
    
    inv <- which(bim_ukb$alt==bim_ref$ref & bim_ukb$ref==bim_ref$alt)
    print(length(inv))
    
    Z[inv] = -Z[inv]
    a <- bim_ukb$alt[inv]
    bim_ukb$alt[inv] <- bim_ukb$ref[inv]
    bim_ukb$ref[inv] <- a
    
    print(table(bim_ukb$alt==bim_ref$alt & bim_ukb$ref==bim_ref$ref))
    
    
    falseid <- which((bim_ukb$alt==bim_ref$alt & bim_ukb$ref==bim_ref$ref)==F)
    
    print(paste("number of falseid:",length(falseid)))
    
    # if(length(falseid)>0){
    #   print("removing false SNPs:")
    #   print(falseid)
    #   Z <- Z[-falseid,]
    #   bim_ref <- bim_ref[-falseid,]
    #   X_ref <- X_ref[,-falseid]
    # }
    # 
    # print(table(bim_ukb$alt==bim_ref$alt & bim_ukb$ref==bim_ref$ref))
    
    LD_ref <- cor(X_ref)
    # library(susieR)
    lambda = estimate_s_rss(Z, LD_ref, n = N.effect)
    print(lambda)
    
    if(regularization == 'regu'){
      LD_ref  = 0.1 * diag(nrow(LD_ref)) + 0.9 * LD_ref
      lambda = estimate_s_rss(Z, LD_ref, n = N.effect)
      print(lambda)
    }
    
    
    ################## 4. ghostknockff
    

    
    start <- Sys.time()
    set.seed(seed)
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD_ref,
      M       = M,
      method  = "sdp" 
    )
    GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    
    end <- Sys.time()
    time4 <- as.numeric(difftime(end, start, units = "secs"))
    cat("ghostknockoff time:", time4, "seconds\n")
    
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej4 <- which(GK.filter$q <= threshold)
      power4[j] = power_cal(rej4, rand)
      fdr4[j] = fdr_cal(rej4, rand)
    }
    
    
    
    ############### 5. GK_simple pseudo sum 
    
    
    start <- Sys.time()
    set.seed(seed)
    GK1_M1_anno_ps  = GK_simple(Z = Z, 
                                R = R, 
                                M = 1, 
                                LD = LD_ref,
                                n = N.effect,
                                ts = 'lasso')
    beta <- GK1_M1_anno_ps$beta
    T_0<-abs(beta[1:p])
    T_k<-abs(matrix(beta[-(1:p)],p,M))
    
    GK.filter<-GhostKnockoff.filter(T_0,T_k)
    
    end <- Sys.time()
    time5 <- as.numeric(difftime(end, start, units = "secs"))
    cat("ghostknockoff time:", time5, "seconds\n")
    
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej5 <- which(GK.filter$q <= threshold)
      power5[j] = power_cal(rej5, rand)
      fdr5[j] = fdr_cal(rej5, rand)
    }
    
    
    
    
    
    ############# 6. GK1 ps with annotation
    
    start <- Sys.time()
    set.seed(seed)
    GK1ps_anno = GK_anno_M(Z, R, M, LD_ref, N.effect)
    
    
    beta <- GK1ps_anno$beta_final
    T_0<-GK1ps_anno$T_0
    T_k<-GK1ps_anno$T_k
    
    GK.filter<-GhostKnockoff.filter(T_0,T_k)
    
    end <- Sys.time()
    time6 <- as.numeric(difftime(end, start, units = "secs"))
    cat("ghostknockoff time:", time6, "seconds\n")
    
    
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej6 <- which(GK.filter$q <= threshold)
      power6[j] = power_cal(rej6, rand)
      fdr6[j] = fdr_cal(rej6, rand)
    }
    
    
    ############ 7. Adaknockoff side with random forest filter
    start <- Sys.time()
    set.seed(seed)
    X_tilde = create.gaussian(X, rep(0,p), LD) # generate knockoff variable
    X_comb = cbind(X,X_tilde)
    
    mdl = cv.glmnet(X_comb,y,alpha=1)
    cvlambda = mdl$lambda.min
    beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
    T0 = abs(beta[1:p])
    T_tilde = abs(beta[(p+1):(2*p)])
    T_max = pmax(T0,T_tilde)
    W1 = T0-T_tilde
    
    # resj = filter_randomForest(W1,z,alpha =alphalist,offset=1)
    resj = filter_randomForest(W1,z,alpha =alphalist,offset=1)
    end <- Sys.time()
    time7 <- as.numeric(difftime(end, start, units = "secs"))
    for(j in 1:len){
      rej7  = resj$rejs[[j]]
      power7[j] = power_cal(rej7, rand)
      fdr7[j] = fdr_cal(rej7, rand)
    }
    
    
    
    ############ 8. Adaknockoff side with EM filter
    # start <- Sys.time()
    # set.seed(seed)
    # X_tilde = create.gaussian(X, rep(0,p), LD) # generate knockoff variable
    # X_comb = cbind(X,X_tilde)
    # 
    # mdl = cv.glmnet(X_comb,y,alpha=1)
    # cvlambda = mdl$lambda.min
    # beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
    # T0 = abs(beta[1:p])
    # T_tilde = abs(beta[(p+1):(2*p)])
    # T_max = pmax(T0,T_tilde)
    # W1 = T0-T_tilde
    # 
    # resj = filter_EM(W1, z, alpha = alphalist, offset = 1, df = 2)
    # end <- Sys.time()
    # time8 <- as.numeric(difftime(end, start, units = "secs"))
    # for(j in 1:len){
    #   rej8  = resj$rejs[[j]]
    #   power8[j] = power_cal(rej8, rand)
    #   fdr8[j] = fdr_cal(rej8, rand)
    # }
    
    
    
    result <- list(region = i,
                   power1 = power1,
                   power2 = power2,
                   power3 = power3,
                   power4 = power4,
                   power5 = power5,
                   power6 = power6,
                   power7 = power7,
                   # power8 = power8,
                   fdr1 = fdr1,
                   fdr2 = fdr2,
                   fdr3 = fdr3,
                   fdr4 = fdr4,
                   fdr5 = fdr5,
                   fdr6 = fdr6,
                   fdr7 = fdr7,
                   # fdr8 = fdr8,
                   time1 = time1,
                   time2 = time2,
                   time3 = time3,
                   time4 = time4,
                   time5 = time5,
                   time6 = time6,
                   time7 = time7,
                   # time8 = time8,
                   X = X,
                   y = y,
                   Xk = X_tilde)
    
    path = paste0(repro_path("true_genotype_legacy_results", "result_external_iter_"), iter, "_n_", N.effect, "_h2e_", h2e, "_alphae_", alphae, "_p1_", p1, "_M_", M, "_chr_", chrid, regularization, ".RData")
    
    save(result, file = path)
    
    
  }
  
}

