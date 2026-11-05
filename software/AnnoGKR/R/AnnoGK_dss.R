# Historical DSS implementation from GK_anno.R; see SOURCE.txt for provenance.
# Function expressions are preserved; comments and formatting are normalized.

#' Annotation-informed GhostKnockoff with SNP-specific sample sizes
#'
#' Reproduce the historical AnnoGK-dss simulation method using summary
#' statistics, an LD matrix, annotations, and SNP-specific sample sizes.
#' @param Z Numeric vector of z-scores, in the same SNP order as LD and R.
#' @param R Numeric annotation matrix with one row per SNP. Nonconstant
#'   columns are standardized internally.
#' @param M Number of knockoff copies, a positive integer.
#' @param LD SNP correlation matrix with dimensions length(Z) by length(Z).
#' @param N.effect Positive effective sample size, either a scalar or a
#'   vector of length length(Z). A scalar is replicated for all SNPs.
#' @param ts Importance-statistic method. Only "lasso" is supported by
#'   the historical DSS helper.
#' @return A list containing original-feature statistics T_0, knockoff
#'   statistics T_k, CV_BEST, lambda_BEST, obj_cor_list, annotation
#'   coefficients lambda_s, and penalty weights Prior.
#' @details This entry point preserves the historical DSS algorithm, including
#'   its tuning and weight transformations. It is not an alias for GK_anno_M.
#'   Set the random seed before fitting. This function does not itself select
#'   SNPs; pass T_0 and T_k to GhostKnockoff.filter and threshold its q values.
#' @seealso [GK_anno_M()], [GhostKnockoff.filter()]
#' @export
GK_anno_dss <- function(Z, R, M = 1, LD, N.effect, ts = "lasso") {
    p = nrow(LD)
    if (length(N.effect) == 1) {
        N.effect = rep(N.effect, p)
    }
    var_para <- apply(R, 2, var)
    N.median = median(N.effect)
    if (all(is.na(R)) | sum(var_para != 0) == 0) {
        print("no informative annotations")
        fit.prelim <- GhostKnockoff.prelim(cor.G = LD, M = M, method = "sdp")
        GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method = "lasso")
        T_0 = GK1_lasso$T_0[[1]]
        T_k = GK1_lasso$T_k[[1]]
        return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA, obj_cor_list = NA, lambda_s = rep(0,
            ncol(R)), Prior = rep(0, 2 * p)))
    }
    R <- R[, var_para != 0]
    R = scale(R)
    fit.prelim <- GhostKnockoff.prelim(cor.G = LD, M = M, method = "sdp")
    GK1_lasso <- GhostKnockoff.lambda(Z, N.effect, fit.prelim, method = ts)
    A = GK1_lasso$A
    modi_para <- GK1_lasso$modi_para
    r_all = GK1_lasso$r_all
    lambda.seq = GK1_lasso$lambda.seq
    lambda = GK1_lasso$lambda
    N.matrix = diag(rep(N.effect, (M + 1))/N.median)
    N.matrix.sqrt <- diag(sqrt(rep(N.effect, M + 1)/N.median))
    A.modi = N.matrix.sqrt %*% A %*% N.matrix.sqrt
    r_all.modi = N.matrix %*% r_all
    temp.A.modi <- (1 - modi_para) * A.modi + diag(modi_para, nrow(A.modi))
    coeff_list <- exp(seq(-2, 2, by = 0.2))
    Prior_final = numeric((M + 1) * p) + 1
    CV_BEST = -Inf
    lambda_BEST = 0
    lambda_s_final = NA
    beta_final = NA
    obj_cor_list <- numeric()
    for (coeff in coeff_list) {
        GKanno_M1_anno_ps <- sol_beta_and_lam_ps_dss(fit.prelim = fit.prelim, Zscore_0 = Z, M = M,
            R = R, temp.A = temp.A.modi, r_all = r_all.modi, N.effect = N.effect, lambda.seq = coeff *
                lambda.seq, lambda_0 = coeff * lambda)
        obj_cor <- cv_func_dss(r_all = r_all, temp.A.modi = temp.A.modi, A.left = fit.prelim$A.left,
            lambda.seq = coeff * lambda.seq, Prior = GKanno_M1_anno_ps$Prior, N.effect = N.effect,
            M = M, nfold = 5)
        obj_cor_list <- append(obj_cor_list, obj_cor)
        if (obj_cor > CV_BEST) {
            lambda_BEST = coeff * lambda
            CV_BEST = obj_cor
            Prior_final = GKanno_M1_anno_ps$Prior
            lambda_s_final = GKanno_M1_anno_ps$lambda_s
            beta_final = GKanno_M1_anno_ps$beta
        }
    }
    T_0 <- abs(beta_final[1:p])
    T_k <- abs(matrix(beta_final[-(1:p)], p, M))
    return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST, obj_cor_list = obj_cor_list,
        lambda_s = lambda_s_final, Prior = Prior_final))
}

#' @noRd
sol_beta_and_lam_ps_dss = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, N.effect, lambda.seq,
    lambda_0, init_scale = NULL, maxiter = 100, verbose = TRUE, init_lam = NULL) {
    if (is.null(init_scale)) {
        init_scale = 0
    }
    n = median(N.effect)
    p = length(Zscore_0)
    Prior = numeric((M + 1) * p) + 1
    fit_basil <- ghostbasil(temp.A, r_all, user.lambdas = lambda.seq * init_scale, alpha = 1, penalty = Prior,
        delta.strong.size = max(1, min(500, length(r_all)/20)), max.strong.size = length(r_all),
        n.threads = 1, use.strong.rule = FALSE)
    beta <- fit_basil$betas[, ncol(fit_basil$betas)]
    beta_old = beta
    lambda_s_old = numeric(ncol(R)) + 100
    lambda_s_all <- numeric()
    objs = c()
    init_lam <- numeric(ncol(R)) + 5
    for (t in 1:100) {
        res <- optim(par = init_lam, fn = obj_fun_gk, method = "BFGS", lambda_0 = lambda_0, beta = beta,
            R = R, n = n, M = M)
        lambda_s = res$par
        print(lambda_s)
        lambda_s_all <- append(lambda_s_all, lambda_s)
        w = weight_standardized(lambda_s, R)
        Prior = rep(w, M + 1)
        fit_basil <- ghostbasil(temp.A, r_all, user.lambdas = lambda.seq * sum(Prior)/((M + 1) *
            p), alpha = 1, penalty = Prior, delta.strong.size = max(1, min(500, length(r_all)/20)),
            max.strong.size = length(r_all), n.threads = 1, use.strong.rule = FALSE)
        beta <- fit_basil$betas[, ncol(fit_basil$betas)]
        obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
        obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
        obj_value = obj_curr - obj_old
        objs = c(objs, obj_curr)
        cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old),
            " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")
        if (mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-04) {
            print(paste0("convergence achieved at ", t, "-th iteration"))
            break
        }
        lambda_s_old = lambda_s
        beta_old = beta
    }
    list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior, lambda_s = lambda_s)
}

#' @noRd
cv_func_dss <- function(r_all, temp.A.modi, A.left, lambda.seq, Prior, N.effect, M = M, nfold = 5) {
    N.median = median(N.effect)
    nA <- N.effect * (nfold - 1)/nfold
    nB <- N.effect/nfold
    temp.left <- sqrt(nB/nA/N.effect) * as.matrix(A.left)
    r_all_A <- r_all + as.matrix(temp.left %*% matrix(rnorm(ncol(temp.left)), ncol(temp.left), 1))
    r_all_B <- (r_all * N.effect - r_all_A * nA)/nB
    N.matrix = diag(rep(N.effect, (M + 1))/N.median)
    r_all_A.modi = N.matrix %*% r_all_A
    fit.basil <- ghostbasil(temp.A.modi, r_all_A.modi, user.lambdas = lambda.seq, alpha = 1, delta.strong.size = max(1,
        min(500, length(r_all_A)/20)), penalty = Prior, max.strong.size = nrow(temp.A.modi), n.threads = 1,
        use.strong.rule = F)
    beta <- fit.basil$betas[, ncol(fit.basil$betas)]
    if (all(beta == 0) == TRUE) {
        return(-Inf)
    }
    Get.f <- function(x) {
        x <- as.matrix(x)
        return(t(x) %*% r_all_B/sqrt(t(x) %*% temp.A.modi %*% x))
    }
    return(Get.f(beta))
}

#' @noRd
GhostKnockoff.lambda <- function(Zscore_0, N.effect, fit.prelim, method = "lasso", type = "fdr",
    M.fwer = 50) {
    Zscore_0 <- as.matrix(Zscore_0)
    N.effect <- as.vector(N.effect)
    N.median <- median(N.effect)
    Zscore_0[is.na(Zscore_0)] <- 0
    N.effect[is.na(N.effect)] <- Inf
    M <- fit.prelim$M
    n.G <- nrow(Zscore_0)
    P.each <- fit.prelim$P.each
    Normal_50Studies <- fit.prelim$Normal_50Studies
    A <- as.matrix(fit.prelim$A)
    A.left <- fit.prelim$A.left
    V.left <- fit.prelim$V.left
    if (type == "fdr") {
        M.rep <- 1
    }
    if (type == "fwer") {
        M.rep <- M.fwer
    }
    for (m in 1:M.rep) {
        Normal_k <- matrix(V.left %*% matrix(rnorm(ncol(V.left)), ncol(V.left), 1), nrow = n.G)
        GK.Zscore_0 <- Zscore_0
        GK.Zscore_k <- as.vector(P.each %*% GK.Zscore_0) + Normal_k
        if (method == "lasso") {
            r <- GK.Zscore_0/sqrt(N.effect)
            r_k <- as.vector(GK.Zscore_k/sqrt(N.effect))
            r_all <- as.matrix(c(r, r_k))
            nfold <- 5
            nA <- N.effect * (nfold - 1)/nfold
            nB <- N.effect/nfold
            temp.left <- sqrt(nB/nA/N.effect) * as.matrix(A.left)
            r_all_A <- r_all + as.matrix(temp.left %*% matrix(rnorm(ncol(temp.left)), ncol(temp.left),
                1))
            r_all_B <- (r_all * N.effect - r_all_A * nA)/nB
            shrink = 0.01
            beta.all <- c()
            parameter.set <- c()
            k <- 1
            N.matrix = diag(rep(N.effect, (M + 1))/N.median)
            N.matrix.sqrt <- diag(sqrt(rep(N.effect, M + 1)/N.median))
            A.modi = N.matrix.sqrt %*% A %*% N.matrix.sqrt
            r_all_A.modi = N.matrix %*% r_all_A
            temp.A.modi <- A.modi + diag(shrink[k], nrow(A))
            fit.basil <- try(ghostbasil(temp.A.modi, r_all_A.modi, alpha = 1, delta.strong.size = max(1,
                min(500, length(r_all_A)/20)), max.strong.size = nrow(temp.A.modi), n.threads = 1,
                use.strong.rule = FALSE), silent = T)
            parameter.set <- rbind(parameter.set, cbind(fit.basil$lmdas, shrink[k]))
            beta.all <- cbind(beta.all, fit.basil$betas)
            Get.f <- function(x) {
                x <- as.matrix(x)
                return(t(x) %*% r_all_B/sqrt(t(x) %*% temp.A.modi %*% x))
            }
            f.lambda <- apply(beta.all, 2, Get.f)
            f.lambda[is.na(f.lambda)] <- -Inf
            parameter <- parameter.set[which.max(f.lambda), ]
            temp.A <- (1 - parameter[2]) * A + diag(parameter[2], nrow(A))
            lambda.all <- fit.basil$lmdas
            lambda <- fit.basil$lmdas[which.max(f.lambda)]
            lambda.seq <- lambda.all[lambda.all > lambda]
            lambda.seq <- c(lambda.seq, lambda)
        }
    }
    return(list(A = A, modi_para = parameter[2], r_all = r_all, lambda.seq = lambda.seq, lambda = lambda))
}
