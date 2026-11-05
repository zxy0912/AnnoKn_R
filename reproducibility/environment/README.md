# Dependencies and software references

The validation environment uses R 4.3.3 and Python 3.11.13 on Linux. [r_dependencies.tsv](r_dependencies.tsv) lists the 43 R packages referenced by the packages, analysis scripts, bundled helpers, and checks, including R's built-in packages and the dynamically loaded `doMC`. It gives checked versions, installation sources, GitHub commits where available, and citation keys. It is a direct-dependency inventory, not a lockfile for all transitive dependencies or a reconstruction of every historical run.

## Installation and checks

For AnnoKn and AnnoGKR alone, follow the [root installation instructions](../../README.md#installation). For all supplied analysis scripts, install the additional packages below from the repository root in R. The explicit versions support R 4.3.3; installing the latest CRAN releases may require a newer R. A compiler toolchain and the system libraries required by the upstream packages are needed for source installations. `doMC` requires a Unix-like system.

```r
# Use your own library; do not modify a shared installation.
lib <- path.expand("~/R/annokn-4.3-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS_USER = lib) # also used by subsequent Rscript checks
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)

deps <- read.delim("reproducibility/environment/r_dependencies.tsv",
                   colClasses = "character", na.strings = "")
# Install Matrix before packages that depend on it.
deps <- deps[order(deps$package != "Matrix"), ]
for (i in which(deps$source == "CRAN" | startsWith(deps$source, "github:"))) {
  p <- deps$package[i]
  installed <- tryCatch(packageVersion(p), error = function(e) NULL)
  if (!is.null(installed) && installed == package_version(deps$version[i])) next
  if (deps$source[i] == "CRAN") {
    remotes::install_version(p, version = deps$version[i], lib = lib,
                             dependencies = NA, upgrade = "never")
  } else {
    remotes::install_github(sub("^github:", "", deps$source[i]), ref = deps$ref[i],
                            lib = lib, dependencies = NA, upgrade = "never")
  }
}
remotes::install_local(".", lib = lib, upgrade = "never")
remotes::install_local("software/AnnoGKR", lib = lib, upgrade = "never")
```

Run these commands from the repository root, using the same R installation and library. If starting a new shell, set `R_LIBS_USER` to the library used above.

```sh
export R_LIBS_USER="$HOME/R/annokn-4.3-library"
Rscript --vanilla reproducibility/environment/check_dependencies.R
Rscript --vanilla reproducibility/tests/test_paths.R reproducibility
Rscript --vanilla reproducibility/tests/test_source_bundle.R reproducibility --load
Rscript --vanilla reproducibility/tests/test_figure2.R reproducibility
```

The dependency checker loads each namespace, checks the minimum versions declared by the two packages, and reports differences from the checked versions. Missing packages, load failures, and violations of declared minimum versions return a nonzero exit status. A different version is reported but is not by itself a failure or evidence of numerical equivalence.

Validation completed with all 43 namespaces loading and all declared minimum versions satisfied. AnnoKn, AnnoGKR, `adaptiveKnockoff`, and `doMC` were installed into a separate library and then loaded; the other dependencies were checked using existing installations. The path/runner, source-bundle loading, and Figure 2 reader checks passed. These checks do not constitute a fresh installation of the entire dependency tree or a rerun of all manuscript experiments. [legacy_package_versions.tsv](legacy_package_versions.tsv) retains the earlier environment inventory, in which `doMC` had not yet been installed.

The supplied Python scripts import only standard-library modules and need no additional pip packages. PLINK v2.00a5LM was used for GTEx and UK Biobank genotype preprocessing. PLINK installation and preprocessing are outside the R dependency check above; obtain PLINK from its [official website](https://www.cog-genomics.org/plink/2.0/).

## Software references

[software_references.bib](software_references.bib) contains BibTeX entries for R, Python, PLINK, the third-party R packages in the dependency inventory, and relevant comparison methods. The inventory's `citation_key` column maps each dependency to an entry. AnnoKn and AnnoGKR are the methods described in the accompanying manuscript. For manuscript citations, prioritize software used for the reported analysis; the utility-package entries are also available for attribution.

| Software / method | Role | Reference / BibTeX key |
| --- | --- | --- |
| R; Python | Analysis and utilities | `R-core`; `Python` |
| PLINK | GTEx and UK Biobank genotype preprocessing | [Chang et al. (2015)](https://doi.org/10.1186/s13742-015-0047-8); `PLINK2-methods`, `PLINK-software` |
| glmnet | Penalized regression | [Friedman et al. (2010)](https://doi.org/10.18637/jss.v033.i01); `R-glmnet` |
| knockoff; SNPknock | Knockoff construction and filtering | `R-knockoff`; [Sesia et al. (2019)](https://doi.org/10.1093/biomet/asy033), `R-SNPknock` |
| adaptiveKnockoff | Adaptive-knockoff comparison | [Ren and Candes (2023)](https://doi.org/10.1214/22-AOAS1663); `R-adaptiveKnockoff` |
| adaptMT | AdaPT comparison | [Lei and Fithian (2018)](https://doi.org/10.1111/rssb.12274); `R-adaptMT` |
| xtune | Annotation-dependent penalty tuning in bundled comparison code | `R-xtune` |
| ghostbasil; GhostKnockoff | Summary-statistic fitting | `R-ghostbasil`; [He et al. (2022)](https://doi.org/10.1038/s41467-022-34932-z), `GhostKnockoff`; [Chen et al. (2024)](https://arxiv.org/abs/2402.12724), `GhostKnockoff-2024` |
| BASIL | Large-scale penalized regression algorithm used by ghostbasil | [Qian et al. (2020)](https://doi.org/10.1371/journal.pgen.1009141); `BASIL` |
| Transfer-knockoff helpers | Transfer-learning comparison | [Li et al. (2022)](https://doi.org/10.1007/s13571-022-00297-y); `Transfer-knockoffs` |
| susieR | Summary-statistic fine-mapping comparison | [Zou et al. (2022)](https://doi.org/10.1371/journal.pgen.1010299); `R-susieR` |
| bigstatsr; genio | Genotype data handling | [Prive et al. (2018)](https://doi.org/10.1093/bioinformatics/bty185), `R-bigstatsr`; `R-genio` |

The remaining numerical, regression, plotting, and data-handling packages are identified individually in the inventory and bibliography. R package references follow the installed packages' `citation()` metadata, with the formal AdaPT publication substituted for its older preprint and the summary-statistic SuSiE reference selected for `susieR`. PLINK has both a software entry with its version and a methods-paper entry, following the [official citation instructions](https://www.cog-genomics.org/plink/2.0/general_usage#cite).
