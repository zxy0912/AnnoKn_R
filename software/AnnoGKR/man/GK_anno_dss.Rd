% Generated from the roxygen comments in R/AnnoGK_dss.R.
\name{GK_anno_dss}
\alias{GK_anno_dss}
\title{Annotation-informed GhostKnockoff with SNP-specific sample sizes}
\usage{
GK_anno_dss(Z, R, M = 1, LD, N.effect, ts = "lasso")
}
\arguments{
\item{Z}{Numeric vector of z-scores, in the same SNP order as LD and R.}
\item{R}{Numeric annotation matrix with one row per SNP. Nonconstant columns
are standardized internally.}
\item{M}{Number of knockoff copies, a positive integer.}
\item{LD}{SNP correlation matrix with dimensions length(Z) by length(Z).}
\item{N.effect}{Positive effective sample size, either a scalar or a vector
of length length(Z). A scalar is replicated for all SNPs.}
\item{ts}{Importance-statistic method. Only \code{"lasso"} is supported by
the historical DSS helper.}
}
\description{
Reproduce the historical AnnoGK-dss simulation method using summary statistics,
an LD matrix, annotations, and SNP-specific sample sizes.
}
\details{
This entry point preserves the historical DSS algorithm, including its tuning
and weight transformations. It is not an alias for \code{GK_anno_M}.
Set the random seed before fitting. This function does not itself select SNPs;
pass \code{T_0} and \code{T_k} to \code{GhostKnockoff.filter} and threshold its
q values.
}
\value{
A list containing original-feature statistics \code{T_0}, knockoff statistics
\code{T_k}, \code{CV_BEST}, \code{lambda_BEST}, \code{obj_cor_list}, annotation
coefficients \code{lambda_s}, and penalty weights \code{Prior}.
}
\seealso{\code{\link{GK_anno_M}}, \code{\link{GhostKnockoff.filter}}}
