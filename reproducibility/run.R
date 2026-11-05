#!/usr/bin/env Rscript
# Usage: Rscript reproducibility/run.R [--check] relative/script.R [script arguments]
runner_args <- base::commandArgs(trailingOnly = TRUE)
check_only <- length(runner_args) && identical(runner_args[[1L]], "--check")
if (check_only) runner_args <- runner_args[-1L]
if (!length(runner_args)) stop("Usage: Rscript reproducibility/run.R [--check] relative/script.R [arguments]", call. = FALSE)
file_arg <- grep("^--file=", base::commandArgs(), value = TRUE)
if (length(file_arg) != 1L) stop("Run this entry point with Rscript.")
runner_file <- sub("^--file=", "", file_arg)
# Rscript may encode spaces in --file as ~+~ when constructing commandArgs().
if (!file.exists(runner_file)) runner_file <- gsub("~+~", " ", runner_file, fixed = TRUE)
repro_dir <- dirname(normalizePath(runner_file, mustWork = TRUE))
options(annokn.repro_dir = repro_dir)
source(file.path(repro_dir, "R", "paths.R"))
config_file <- Sys.getenv("ANNOKN_REPRO_CONFIG", "")
if (!nzchar(config_file) && file.exists(file.path(repro_dir, "config.local.R")))
  config_file <- file.path(repro_dir, "config.local.R")
config <- repro_load_config(repro_dir, config_file)
script <- normalizePath(file.path(repro_dir, runner_args[[1L]]), mustWork = TRUE)
if (!startsWith(script, paste0(repro_dir, "/")) ||
    !grepl("/(simulations|real_data)/.*\\.(R|py)$", script))
  stop("Select an analysis script under reproducibility/simulations or real_data.")
check <- repro_preflight(config, script)
if (!check$ok) stop("Missing configured inputs/sources. Edit config.local.R or set ANNOKN_REPRO_CONFIG. No analysis was run.", call. = FALSE)
if (check_only) {
  cat("Path check passed. No analysis was run and no output directories were created.\n")
} else {
  for (key in check$output_keys) {
    path <- config$paths[[key]]
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE))
      stop("Could not create output directory: ", path)
  }
  repro_path <- repro_path_function(config)
  for (ref in repro_script_paths(script)) {
    if (!ref$key %in% check$output_keys) next
    parent <- dirname(repro_path(ref$key, ref$relative))
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE))
      stop("Could not create output subdirectory: ", parent)
  }
  analysis_args <- runner_args[-1L]
  # Historical drivers expect only their own arguments (e.g. replicate/chromosome).
  commandArgs <- local({
    supplied <- analysis_args
    all_args <- c(base::commandArgs()[1L], paste0("--file=", script), "--args", supplied)
    function(trailingOnly = FALSE) if (trailingOnly) supplied else all_args
  })
  setwd(dirname(repro_dir))
  if (grepl("\\.py$", script)) {
    resolved_file <- tempfile("annokn-resolved-paths-", fileext = ".tsv")
    write.table(data.frame(key = names(config$paths), path = unname(config$paths)),
                resolved_file, sep = "\t", row.names = FALSE, quote = TRUE)
    status <- tryCatch(system2(Sys.getenv("ANNOKN_PYTHON", "python3"),
      shQuote(c(file.path(repro_dir, "python", "run_analysis.py"),
                "--paths", resolved_file, script, analysis_args))),
      finally = unlink(resolved_file))
    if (status != 0L) quit(status = status)
  } else {
    sys.source(script, envir = .GlobalEnv)
  }
}
