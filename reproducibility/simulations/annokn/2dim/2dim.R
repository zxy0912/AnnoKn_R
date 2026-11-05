# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)

args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

i <- args[1]

####################################

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

rho = 0.5

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)

cortype = 1

print(paste0("######################################### iteration ", i, "  ######################################### "))

seed = 1000 * i

n = 1000
p = 1600 # 40-by-40
k = 200
dist_para = 3

alphalist <- seq(0.3,0.05,-0.05)
Sigma = diag(rep(1,p))
amp = 25
type = 'in'
binary = 'binary'

for(cortype in c(3)){
  for (amp in c(25)){
    #for(type in c(0.5, 1, 2, 'in')){
    for(type in c(1, 0.5, 2, 'in')){
      for(binary in c('','binary')){
        
        set.seed(seed)
        
        power1 = numeric(len)
        power2 = numeric(len)
        power3 = numeric(len)
        power4 = numeric(len)
        power5 = numeric(len)
        power6 = numeric(len)
        power7 = numeric(len)
        fdr1 = numeric(len)
        fdr2 = numeric(len)
        fdr3 = numeric(len)
        fdr4 = numeric(len)
        fdr5 = numeric(len)
        fdr6 = numeric(len)
        fdr7 = numeric(len)
        
        # cv_matrix <- numeric(61)
        # cv_matrix_em <- numeric(61)
        
        all_res <- data.frame()
        
        ####################################
        ## determining the signals
        ####################################
        nonzero = c()
        
        cor_hypothesis <- expand.grid(x = 1:40, y = 1:40)
        cor_hypothesis$index <- 1:nrow(cor_hypothesis)
        # causal_region <- subset(cor_hypothesis, x <= 20 & y <= 20)
        
        if(type == 1){
          causal_region <- subset(cor_hypothesis, x + y <= 30)
        }else if(type == 2){
          causal_region <- subset(cor_hypothesis, x^2 + y^2 <= 600)
        }else if(type == 0.5){
          causal_region <- subset(cor_hypothesis, sqrt(x) + sqrt(y) <= sqrt(50))
        }else if(type == 'in'){
          causal_region <- subset(cor_hypothesis, x <= 20 & y <= 20)
        }
        
        causal_region$prob <- 1 / (causal_region$x + causal_region$y)
        #causal_region$prob <- 1/((causal_region$x)^d + (causal_region$y)^d)
        causal_region$prob <- causal_region$prob / sum(causal_region$prob)  
        selected_rows <- sample(1:nrow(causal_region), size = k, prob = causal_region$prob)
        nonzero <- causal_region$index[selected_rows[order(selected_rows)]]
        
        
        
        # alphalist = seq(0.3,0.01,-0.01)
        # covariance matrix
        # type1: 
        if(cortype == 1){
          Sigma = diag(p)
        }else{
          Sigma = exp(-dist_para*pdist(cor_hypothesis[, c(1,2)], metric = "euclidean", p = 2)^2)
        }
        
        
        # off_diag_values <- Sigma[row(Sigma) != col(Sigma)]
        # summary(off_diag_values)
        k = length(nonzero)
        beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
        y.sample = function(X) rbinom(n,1,exp(X %*% beta0)/(1+exp(X %*% beta0)))
        
        ####################################
        ## Generating data
        ####################################
        X = matrix(rnorm(n*p),n) %*% chol(Sigma)
        y = y.sample(X)
        
        # standardize the matrix
        # X = scale(X)
        # Xk = scale(Xk)
        # y = y / sd(y)
        
        
        ############################ making the plot
        
        indi <- numeric(p)
        indi[nonzero] = 1
        col_vec <- ifelse(indi == 1, "black", "grey")
        
        plot1 <- plot(x = cor_hypothesis$x,
                      y = cor_hypothesis$y,
                      col = col_vec,
                      pch = 16,
                      xlab = "r",
                      ylab = "s",
                      main = "Scenario 1")
        
        print(plot1)
        
        
        z = as.matrix(cor_hypothesis[,1:2])
        
        if(binary == 'binary'){
          z[, 1] = ifelse(z[,1] <= max(causal_region$y), 0 , 1)
          z[, 2] = ifelse(z[,2] <= max(causal_region$y), 0 , 1)
        }
        
        R_inf <- scale(as.matrix(z))
        p         <- nrow(R_inf)
        R <- as.matrix(R_inf)
        
        
        
        
        
        ############# 1. Vanilla  knockoff
        
        start = Sys.time()
        Xk = create.gaussian(X,rep(0,p),Sigma)
        end_check1 = Sys.time()
        mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
        cvlambda = mdl$lambda.min
        beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
        T = abs(beta[1:p])
        T_tilde = abs(beta[(p+1):(2*p)])
        T_max = pmax(T,T_tilde)
        W1 = T-T_tilde
        end_check2 = Sys.time()
        
        time1 <- as.numeric(difftime(end_check2, start, units = "secs"))
        time_check1 <- as.numeric(difftime(end_check1, start, units = "secs"))
        
        for (j in 1:len) {
          alpha = alphalist[j]
          tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
          rej1 = as.numeric(which(W1>=tau))
          power1[j] = power_cal(rej1, nonzero)
          fdr1[j] = fdr_cal(rej1, nonzero)
          
        }
        
        
        ############ 2. Adaknockoff with glm
        
        start = Sys.time()
        resj = filter_glm(W1,z,alpha =alphalist,offset=1)
        end = Sys.time()
        time2 <- as.numeric(difftime(end, start, units = "secs")) + time1
        
        for(j in 1:len){
          rej2  = resj$rejs[[j]]
          power2[j] = power_cal(rej2, nonzero)
          fdr2[j] = fdr_cal(rej2, nonzero)
        }
        
        
        ########### 3. Adaknockoff with gam
        time3 = 0
        # if(binary == ''){
        #   start = Sys.time()
        #   resj = filter_gam(W1,z,alpha =alphalist,offset=1)
        #   end = Sys.time()
        #   time3 <- as.numeric(difftime(end, start, units = "secs")) + time1
        #   
        #   for(j in 1:len){
        #     rej3  = resj$rejs[[j]]
        #     power3[j] = power_cal(rej3, nonzero)
        #     fdr3[j] = fdr_cal(rej3, nonzero)
        #   }
        # }
        
        ######### 4. Adaknockoff with random forest
        start = Sys.time()
        resj = filter_randomForest(W1,z,alpha =alphalist,offset=1)
        end = Sys.time()
        time4 <- as.numeric(difftime(end, start, units = "secs")) + time1
        
        for(j in 1:len){
          rej4  = resj$rejs[[j]]
          power4[j] = power_cal(rej4, nonzero)
          fdr4[j] = fdr_cal(rej4, nonzero)
        }
        
        ######### 5. Adaknockoff with EM algorithm
        
        
        if(binary == ''){
          start = Sys.time()
          resj =filter_EM(W1,z,alpha =alphalist,offset=1,cutoff = 0)
          end = Sys.time()
          time5 <- as.numeric(difftime(end, start, units = "secs")) + time1
          for(j in 1:len){
            rej5  = resj$rejs[[j]]
            power5[j] = power_cal(rej5, nonzero)
            fdr5[j] = fdr_cal(rej5, nonzero)
          }
        }else{
          start = Sys.time()
          resj =filter_EM(W1,z,alpha =alphalist,offset=1,cutoff = 0)
          end = Sys.time()
          time5 <- as.numeric(difftime(end, start, units = "secs")) + time1
          for(j in 1:len){
            rej5  = resj$rejs[[j]]
            power5[j] = power_cal(rej5, nonzero)
            fdr5[j] = fdr_cal(rej5, nonzero)
          }
        }
        
        
        
        
        
        ######### 6. Knockoff-anno
        set.seed(seed)
        R <- scale(as.matrix(z))
        start = Sys.time()
        result_annokn <- knockoff_anno_improved(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
        end = Sys.time()
        time6 <- as.numeric(difftime(end, start, units = "secs")) + time_check1
        beta = result_annokn$beta
        
        lambda_vals <- unlist(result_annokn$lambda_s)
        
        T = abs(beta[1:p])
        T_tilde = abs(beta[(p+1):(2*p)])
        T_max = pmax(T,T_tilde)
        W3 = T-T_tilde
        
        for (j in 1:len) {
          
          alpha = alphalist[j]
          tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
          rej6 = as.numeric(which(W3>=tau))
          power6[j] = power_cal(rej6, nonzero)
          fdr6[j] = fdr_cal(rej6, nonzero)
          
        }
        
        
        ############# 7. knockoff_simple
        
        set.seed(seed)
        start <- Sys.time()
        result_annokn_lite = knockoff_simple(X = X, Xk = Xk, y = y, R = R)
        end <- Sys.time()
        time7 <- as.numeric(difftime(end, start, units = "secs")) + time_check1
        beta = result_annokn_lite$beta
        
        
        T = abs(beta[1:p])
        T_tilde = abs(beta[(p+1):(2*p)])
        T_max = pmax(T,T_tilde)
        W2 = T-T_tilde
        
        for (j in 1:len) {
          
          alpha = alphalist[j]
          tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
          rej7 = as.numeric(which(W2>=tau))
          power7[j] = power_cal(rej7, nonzero)
          fdr7[j] = fdr_cal(rej7, nonzero)
          
        }
        
        result = list(power1 = power1, # knockoff
                      power2 = power2, # AdaKn (GLM)
                      power3 = power3, # AdaKn (GAM) 
                      power4 = power4, # AdaKn (RF) 
                      power5 = power5, # AdaKn (EM)
                      power6 = power6, # knockoff anno
                      power7 = power7, # knockoff simple
                      lambda_vals = lambda_vals,
                      fdr1 = fdr1,
                      fdr2 = fdr2,
                      fdr3 = fdr3,
                      fdr4 = fdr4,
                      fdr5 = fdr5,
                      fdr6 = fdr6,
                      fdr7 = fdr7,
                      time1 = time1,
                      time2 = time2,
                      time3 = time3,
                      time4 = time4,
                      time5 = time5,
                      time6 = time6,
                      time7 = time7)
        
        
        path = paste0(repro_path("annokn_2d_results", "result_"), binary, "_amp_",amp, "_type_", type, "_cor_", cortype, "_", i,".RData")
        
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
    }
  }
}
