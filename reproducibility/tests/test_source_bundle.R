# Check all active analysis source paths using repository defaults (no local profile).
# Add --load to execute source files; this requires their R package dependencies.
args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript tests/test_source_bundle.R reproducibility [--load]")
base <- normalizePath(args[1], mustWork = TRUE)
source(file.path(base, "R", "paths.R"))
config <- repro_load_config(base)
path <- repro_path_function(config)
keys <- c("annokn_source", "annogk_source", "comparison_utils", "revision_helpers")
stopifnot(all(config$spec$root[match(keys, config$spec$key)] == "repository"))
scripts <- unlist(lapply(c("simulations", "real_data"), function(folder)
  list.files(file.path(base, folder), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)))
sources <- character()
for (script in scripts) {
  for (ref in repro_script_paths(script)) {
    if (!ref$key %in% keys) next
    file <- path(ref$key, ref$relative)
    if (!grepl("[.]R$", file)) file <- paste0(file, ".R")
    stopifnot(startsWith(file, paste0(base, "/")), file.exists(file))
    sources <- c(sources, file)
  }
}
sources <- sort(unique(sources))
stopifnot(length(sources) == 16L) # 11 historical files and 5 revision helpers
# These source bundles must not hide another external source dependency.
inspect <- function(x) {
  if (missing(x)) return(invisible(NULL))
  if (is.call(x) && (identical(x[[1]], as.name("source")) ||
                     identical(x[[1]], as.name("sys.source"))))
    stop("Nested source dependency found; add it to the bundle and audit: ", deparse(x)[1])
  if (is.call(x) || is.expression(x) || is.pairlist(x))
    for (child in as.list(x)) inspect(child)
}
for (file in sources) inspect(parse(file, keep.source = FALSE))
cat("PASS:", length(sources), "active source files resolve inside the checkout; no nested source calls.\n")

if ("--load" %in% args) {
  # Same source order as the Figure 1 driver; separate environment from revision helpers.
  legacy <- new.env(parent = .GlobalEnv)
  utilities <- strsplit(config$spec$expected_content[config$spec$key == "comparison_utils"], "; ")[[1]]
  for (file in utilities) source(path("comparison_utils", file), local = legacy)
  source(path("annokn_source", "KF_anno.R"), local = legacy)
  expected <- c("knockoff_anno_improved", "knockoff_simple", "power_cal", "fdr_cal",
                "filter_EM", "filter_gam", "filter_glm", "filter_randomForest",
                "sampleHMM", "knockoffHMM", "adapt_gam", "cv.glmnet", "knockoff.threshold")
  stopifnot(all(vapply(expected, function(n) exists(n, envir = legacy, mode = "function"), logical(1))))
  cat("PASS: Figure 1 sources load and required model/comparator functions resolve.\n")
  source(path("annogk_source", "GK_anno.R"), local = legacy)
  source(path("annogk_source", "GhostKnockoff.R"), local = legacy)
  stopifnot(exists("GK_anno_dss", envir = legacy, mode = "function"),
            exists("GhostKnockoff.prelim", envir = legacy, mode = "function"))
  revision <- new.env(parent = .GlobalEnv)
  for (file in c("filter_EM.R", "filter_gam.R", "filter_glm.R",
                 "filter_randomForest.R", "transfer_functions.R"))
    source(path("revision_helpers", file), local = revision)
  cat("PASS: historical AnnoGK/GhostKnockoff and separate revision helper sources load.\n")
}
