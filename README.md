# AnnoKn and AnnoGKR

This repository provides two independently installable R packages for annotation-informed variable selection.

| Package | Location | Purpose | Main functions |
| --- | --- | --- | --- |
| **AnnoKn** | Repository root | Individual-level data with functional annotations | `AnnoKn()`, `AnnoKn_lite()` |
| **AnnoGKR** | [`software/AnnoGKR/`](software/AnnoGKR/) | Summary statistics with an LD matrix and annotations | `GK_anno_M()`, `GK_simple()`, `GK_anno_dss()` |

Manuscript simulations, real-data analyses, and their figure/table correspondence are described in [`reproducibility/README.md`](reproducibility/README.md). Data sources and access instructions are listed in its [Data availability section](reproducibility/README.md#data-availability).

## Installation

Install from R. AnnoGKR requires `ghostbasil`, which must be installed separately from GitHub.

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Individual-level package
remotes::install_github("shi-ym19/AnnoKn_R", ref = "jasa-reproducibility")

# Summary-statistic package
remotes::install_github("szcf-weiya/ghostbasil", upgrade = "never")
remotes::install_github(
  "shi-ym19/AnnoKn_R", subdir = "software/AnnoGKR",
  ref = "jasa-reproducibility", upgrade = "never"
)
```

For a local checkout, use `remotes::install_local(".")` and `remotes::install_local("software/AnnoGKR")` after installing dependencies. Package dependencies are declared in the two `DESCRIPTION` files. A Linux Conda environment specification is available in [`software/AnnoGKR/AnnoGK.yml`](software/AnnoGKR/AnnoGK.yml).

The [dependency guide](reproducibility/environment/README.md) provides checked package versions, installation commands for the analysis scripts, a dependency checker, and software references.

## Usage

**AnnoKn:** supply an individual-level design matrix `X`, its knockoff matrix `Xk`, response `y`, and annotation matrix `R` (one row per feature).

```r
set.seed(1000)
fit <- AnnoKn::AnnoKn(X = X, Xk = Xk, y = y, R = R, attempts = c(0))
p <- ncol(X)
W <- abs(fit$beta[seq_len(p)]) - abs(fit$beta[p + seq_len(p)])
selected <- which(W >= knockoff::knockoff.threshold(W, fdr = 0.1, offset = 1))

# AnnoKn-lite uses the same input matrices
fit_lite <- AnnoKn::AnnoKn_lite(X = X, Xk = Xk, y = y, R = R)
```

**AnnoGKR:** supply SNP-level statistics `Z`, LD matrix `LD`, annotations `R`, and effective sample size `N.effect`. SNPs must have the same order across inputs. The following uses AnnoGK-dss:

```r
set.seed(1000)
fit <- AnnoGKR::GK_anno_dss(Z = Z, R = R, M = 1, LD = LD, N.effect = N.effect)
q <- AnnoGKR::GhostKnockoff.filter(fit$T_0, fit$T_k)$q
selected <- which(q <= 0.1)
```

`GK_anno_M()` and `GK_simple()` provide the other annotation-informed variants. For direct baseline calls, use `GhostKnockoff.prelim(..., method = "sdp")` and `GhostKnockoff.fit(..., method = "lasso")`. Use package-qualified function names when loading both packages, as some exported names overlap. Function arguments and return values are documented in the installed R help pages.

## Licenses

AnnoKn uses the [MIT license](LICENSE.md). AnnoGKR retains its GPL-3 license; upstream attribution and package modifications are recorded in [`software/AnnoGKR/SOURCE.txt`](software/AnnoGKR/SOURCE.txt).
