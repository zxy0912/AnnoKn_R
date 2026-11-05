# Base-R integration checks; no statistical models, external data, or packages.
args <- commandArgs(TRUE)
base <- normalizePath(args[1], mustWork = TRUE)
source(file.path(base, "R", "paths.R"))
scratch <- tempfile("annokn path tests ")
dir.create(scratch)
relocated <- file.path(scratch, "moved checkout", "reproducibility")
dir.create(file.path(relocated, "R"), recursive = TRUE)
dir.create(file.path(relocated, "config"))
dir.create(file.path(relocated, "simulations"))
dir.create(file.path(relocated, "python"))
file.copy(file.path(base, "python", "run_analysis.py"), file.path(relocated, "python"))
file.copy(file.path(base, "run.R"), relocated)
file.copy(file.path(base, "R", "paths.R"), file.path(relocated, "R"))
file.copy(file.path(base, "config", "paths.tsv"), file.path(relocated, "config"))
input <- file.path(scratch, "input files")
output <- file.path(scratch, "output files")
dir.create(input)
writeLines("test only", file.path(input, "input.txt"))
dput(list(output_root = output, paths = list(gtex_expression = input)),
     file.path(relocated, "config.local.R"))
script <- file.path(relocated, "simulations", "fixture.R")
writeLines(c(
  'stopifnot(identical(commandArgs(TRUE), c("7", "argument with spaces")))',
  'stopifnot(grepl("fixture.R$", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))))',
  'stopifnot(identical(readLines(repro_path("gtex_expression", "input.txt")), "test only"))',
  'saveRDS(commandArgs(TRUE), repro_path("revision_plots", "nested/probe.rds"))'
), script)
old_env <- Sys.getenv("ANNOKN_REPRO_CONFIG", unset = NA_character_)
Sys.unsetenv("ANNOKN_REPRO_CONFIG")
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(check = FALSE) {
  opts <- c("--vanilla", file.path(relocated, "run.R"),
            if (check) "--check", "simulations/fixture.R", "7", "argument with spaces")
  log <- suppressWarnings(system2(rscript, shQuote(opts), stdout = TRUE, stderr = TRUE))
  code <- attr(log, "status")
  if (is.null(code)) code <- 0L
  list(status = code, log = log)
}
checked <- run(TRUE)
if (checked$status != 0L) cat(checked$log, sep = "\n")
stopifnot(checked$status == 0L, !dir.exists(output))
executed <- run()
if (executed$status != 0L) cat(executed$log, sep = "\n")
stopifnot(executed$status == 0L,
          identical(readRDS(file.path(output, "plots/revision/nested/probe.rds")),
                    c("7", "argument with spaces")))
cfg <- repro_load_config(relocated, file.path(relocated, "config.local.R"))
path <- repro_path_function(cfg)
stopifnot(identical(path("gtex_expression", ""), paste0(input, "/")),
          identical(path("gtex_expression", "chr_"), paste0(input, "/chr_")))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
must_fail(path("misspelled_key"))
must_fail(path("plots", "../escape"))
bad <- file.path(scratch, "bad.R")
dput(list(paths = list(misspelled_key = input)), bad)
must_fail(repro_load_config(relocated, bad))
# Actual missing data and missing source files fail before executing the fixture.
unlink(input, recursive = TRUE)
unlink(output, recursive = TRUE)
missing <- run(TRUE)
stopifnot(missing$status != 0L, !dir.exists(output))
writeLines('source(repro_path("annokn_source", "KF_anno.R"))', script)
missing_source <- run(TRUE)
stopifnot(missing_source$status != 0L,
          any(grepl("MISSING SOURCE", missing_source$log)), !dir.exists(output))
# Python shares the R configuration, preserves script argv/__file__, and checks
# missing roots before execution. Relocated paths contain spaces.
dir.create(input)
writeLines("test only", file.path(input, "input.txt"))
py <- file.path(relocated, "simulations", "fixture.py")
writeLines(c(
  'from pathlib import Path', 'import sys',
  'assert sys.argv[1:] == ["7", "argument with spaces"]',
  'assert Path(__file__).name == "fixture.py"',
  'assert Path(repro_path("gtex_expression", "input.txt")).read_text().strip() == "test only"',
  'Path(repro_path("revision_plots", "nested/python.txt")).write_text("passed")'
), py)
run_python <- function(check = FALSE) {
  opts <- c("--vanilla", file.path(relocated, "run.R"), if (check) "--check",
            "simulations/fixture.py", "7", "argument with spaces")
  log <- suppressWarnings(system2(rscript, shQuote(opts), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (!is.null(status)) cat(log, sep = "\n")
  if (is.null(status)) 0L else status
}
stopifnot(run_python(TRUE) == 0L, !dir.exists(output), run_python() == 0L,
          identical(readLines(file.path(output, "plots/revision/nested/python.txt"), warn=FALSE), "passed"))
unlink(input, recursive = TRUE)
unlink(output, recursive = TRUE)
stopifnot(run_python(TRUE) != 0L, !dir.exists(output))
if (is.na(old_env)) Sys.unsetenv("ANNOKN_REPRO_CONFIG") else Sys.setenv(ANNOKN_REPRO_CONFIG = old_env)
unlink(scratch, recursive = TRUE)
cat("PASS: R/Python relocation, spaces, argument forwarding, script identity, output subdirectories, dry-run behavior, unknown keys, and missing data/source checks.\n")
