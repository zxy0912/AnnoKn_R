if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
args = commandArgs(TRUE)
arg = as.numeric(args[1])


library(adaptiveKnockoff)
library(knockoff)
library(gam)
library(glmnet)
library(AnnoKn)

source(repro_path("revision_helpers", "transfer_functions.R"))
type = 'single'
error = 'corrected'

# Problem parameters
num_env = 3 ##the first env is the target env
overlap_all = c(0.2, 0.8)
ns = rep(800, num_env)

p = 500
q = 0.1
rhos = rep(0.5, num_env)
amplitude = 3.5

N = 1

k = 60
k_overlap = floor(k*overlap_all)

beta = matrix(rep(0,num_env*p), nrow = p)
beta[1:k,1] = 1
for (i in 1:(num_env - 1)) {
  row_indices = 1:k + k - k_overlap[i]
  beta[row_indices, i + 1] = 1
}
beta = beta* amplitude / sqrt(ns[1])
set.seed(1) ## make sure our model is not changing
permu = sample(1:p)
beta = beta[permu,]
nonnulls = which(beta[,1] != 0)


fdps = NULL
powers = NULL
fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)



set.seed(arg)
method_one = c("pool", "vanilla", "adaptive", "gam", paste("linear_com", 0.1*(1:5), sep = "_"), "weighted_lasso", "AnnoKn")
methods = NULL


for (iter in 1:N){

  Xs = list()
  ys = list()
  for (i in 1:num_env){
    Sigma = toeplitz(rhos[i]^(0:(p-1)))
    X = matrix(rnorm(ns[i]*p),ns[i]) %*% chol(Sigma)
    Xs[[i]] = X
    y = X %*% beta[,i] + rnorm(ns[i])
    ys[[i]] = y
  }

  nonnulls = which(beta[,1]>0)

  ## create knockoffs
  X_ks = list()
  for (i in 1:num_env){
    X_ks[[i]] = create.gaussian(Xs[[i]], rep(0,p), toeplitz(rhos[i]^(0:(p-1))))
  }

  ## stack
  X_stack = NULL
  X_k_stack = NULL
  y_stack = NULL
  for (i in 1:num_env){
    X_stack = rbind(X_stack, Xs[[i]])
    X_k_stack = rbind(X_k_stack, X_ks[[i]])
    y_stack = c(y_stack,ys[[i]])
  }
  W_stack = stat.lasso_coefdiff(X_stack, X_k_stack, y_stack)
  thres_stack = knockoff.threshold(W_stack, fdr = q)
  selected_stack = which(W_stack >= thres_stack)
  fdp_stack = fdp(selected_stack, nonnulls)
  power_stack = power(selected_stack,nonnulls)

  ## Vanilla knockoffs
  W_vanilla = stat.lasso_coefdiff(Xs[[1]], X_ks[[1]], ys[[1]])
  thres_vanilla = knockoff.threshold(W_vanilla, fdr = q)
  selected_vanilla = which(W_vanilla >= thres_vanilla)
  fdp_vanilla = fdp(selected_vanilla, nonnulls)
  power_vanilla = power(selected_vanilla,nonnulls)

  ## Our method -- transfer learning -- adaptive knockoffs
  fdp_side = NA
  power_side = NA
  X_stack_other = NULL
  X_stack_k_other = NULL
  y_stack_other = NULL
  for (k2 in 2:num_env){
    X_stack_other = rbind(X_stack_other, Xs[[k2]])
    X_stack_k_other = rbind(X_stack_k_other, X_ks[[k2]])
    y_stack_other = c(y_stack_other, ys[[k2]])
  }
  if(error == 'corrected'){
    W_abs_other = abs(stat.lasso_coefdiff(X_stack_other, X_stack_k_other, y_stack_other))
  }else{
    W_abs_other = abs(stat.lasso_coefdiff(X_stack, X_k_stack, y_stack))
  }

  tryCatch({
    adap = filter_EM(W_vanilla, W_abs_other, alpha = q)
    selected_side = adap$rejs[[1]]
    fdp_side = fdp(selected_side, nonnulls)
    power_side = power(selected_side, nonnulls)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  ## Our method -- transfer learning -- adaptive knockoffs -- gam
  fdp_gam = NA
  power_gam = NA
  adap = filter_gam(W_vanilla, W_abs_other, alpha = q)
  selected_gam = adap$rejs[[1]]
  fdp_gam = fdp(selected_gam, nonnulls)
  power_gam = power(selected_gam, nonnulls)

  ## Our method -- transfer learning -- linear_combination
  fdp_trans2s = NULL
  power_trans2s = NULL

  for (theta in 0.1*(1:5)){
    W_trans2 = compute_transfer_stats_linear_combination(Xs, X_ks, ys, theta = theta)
    thres_trans2 = knockoff.threshold(W_trans2, fdr = q)
    selected_trans2 = which(W_trans2 >= thres_trans2)
    fdp_trans2 = fdp(selected_trans2, nonnulls)
    power_trans2 = power(selected_trans2,nonnulls)
    fdp_trans2s = c(fdp_trans2s, fdp_trans2)
    power_trans2s = c(power_trans2s, power_trans2)
  }

  ## Our method -- transfer learning -- weighted lasso
  W_trans3 = compute_transfer_stats_weighted_lasso(Xs, X_ks, ys)
  thres_trans3 = knockoff.threshold(W_trans3, fdr = q)
  selected_trans3 = which(W_trans3 >= thres_trans3)
  fdp_trans3 = fdp(selected_trans3, nonnulls)
  power_trans3 = power(selected_trans3, nonnulls)


  ## AnnoKn
  if(type == 'all'){
    W_annotation_matrix = matrix(NA, nrow = p, ncol = num_env - 1)
    for (i in 2:num_env){
      W_env_i = stat.lasso_coefdiff(Xs[[i]], X_ks[[i]], ys[[i]])
      W_annotation_matrix[, i - 1] = abs(W_env_i)
    }
  }else{
    W_annotation_matrix = matrix(W_abs_other, ncol = 1)
  }

  W_annotation_matrix = scale(W_annotation_matrix)

  result_annokn <- AnnoKn(X = Xs[[1]], Xk = X_ks[[1]], y = ys[[1]],
                          attempts = c(0), R = W_annotation_matrix)
  lambda_s = result_annokn$lambda_s
  beta_annokn = result_annokn$beta
  T = abs(beta_annokn[1:p])
  T_tilde = abs(beta_annokn[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W_annokn = T-T_tilde

  thres_annokn = knockoff.threshold(W_annokn,fdr = q,offset = 1)
  selected_annokn = as.numeric(which(W_annokn >= thres_annokn))

  fdp_annokn = fdp(selected_annokn, nonnulls)
  power_annokn = power(selected_annokn, nonnulls)

  fdps = c(fdps, c(fdp_stack, fdp_vanilla, fdp_side, fdp_gam, fdp_trans2s, fdp_trans3, fdp_annokn))
  powers = c(powers, c(power_stack, power_vanilla, power_side, power_gam, power_trans2s, power_trans3, power_annokn))
  methods = c(methods, method_one)

}

save(arg, fdps,powers, methods, lambda_s, file = paste(repro_path("transfer_results", "new_output_"),type, error, arg,".Rdata", sep = ""))
