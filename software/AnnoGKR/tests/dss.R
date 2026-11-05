# Exercise the installed namespace: no source files or external test packages.
library(AnnoGKR)
stopifnot("GK_anno_dss" %in% getNamespaceExports("AnnoGKR"))
stopifnot(!any(c("sol_beta_and_lam_ps_dss", "cv_func_dss", "GhostKnockoff.lambda") %in%
              getNamespaceExports("AnnoGKR")))

p <- 24L
LD <- toeplitz(0.5^(0:(p - 1L)))
set.seed(910)
Z <- as.vector(t(chol(LD)) %*% rnorm(p)) +
  as.vector(LD[, 1:6] %*% rep(6, 6))
R <- cbind(position = seq_len(p), binary = rep(c(0, 1), length.out = p))

fit_dss <- function(M, N.effect) {
  set.seed(1000)
  printed <- capture.output(fit <- AnnoGKR::GK_anno_dss(
    Z = Z, R = R, M = M, LD = LD, N.effect = N.effect))
  stopifnot(length(fit$T_0) == p, identical(dim(fit$T_k), c(p, as.integer(M))),
            length(fit$Prior) == (M + 1L) * p,
            length(fit$lambda_s) == ncol(R),
            all(is.finite(c(fit$T_0, fit$T_k, fit$lambda_s, fit$Prior))),
            all(fit$Prior > 0), is.finite(fit$CV_BEST),
            fit$lambda_BEST > 0, any(fit$T_0 > 0))
  q <- AnnoGKR::GhostKnockoff.filter(fit$T_0, fit$T_k)$q
  stopifnot(length(q) == p, all(is.finite(q)), all(q >= 0 & q <= 1))
  fit
}

# A scalar sample size must be equivalent to its per-SNP expansion.
scalar <- fit_dss(1L, 5000)
expanded <- fit_dss(1L, rep(5000, p))
stopifnot(identical(scalar, expanded))

# Multiple copies with genuinely varying per-SNP sample sizes.
varying <- fit_dss(5L, seq(3000, 6000, length.out = p))
stopifnot(identical(varying, fit_dss(5L, seq(3000, 6000, length.out = p))))
