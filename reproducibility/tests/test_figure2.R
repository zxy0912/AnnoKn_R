# Integration checks for Figure 2's complete-input requirement and exported means.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: Rscript tests/test_figure2.R reproducibility")
base <- normalizePath(args[1], mustWork = TRUE)
work <- tempfile("figure2 reader tests ")
dir.create(work)
fixture <- file.path(work, "checkout", "reproducibility")
for (folder in c("R", "config", "simulations/annogk"))
  dir.create(file.path(fixture, folder), recursive = TRUE, showWarnings = FALSE)
for (file in c("run.R", "R/paths.R", "config/paths.tsv", "simulations/annogk/figure2_result.R"))
  stopifnot(file.copy(file.path(base, file), file.path(fixture, file)))
data_dir <- file.path(work, "result files")
dir.create(data_dir)
output_dir <- file.path(work, "plot outputs")
config_file <- file.path(work, "profile.R")
dput(list(output_root = output_dir, paths = list(annogk_results = data_dir)), config_file)
Sys.setenv(ANNOKN_REPRO_CONFIG = config_file)
alpha <- seq(0.4, 0.1, -0.05)
# Vary replicate, method, threshold and heritability; replicate mean is known (50.5).
for (h in c(0.05, 0.1, 0.2)) for (i in 1:100) {
  result <- list()
  for (j in 1:11) {
    result[[paste0("power", j)]] <- i / 1000 + j / 100 + alpha / 10 + h
    result[[paste0("fdr", j)]] <- i / 2000 + j / 200 + alpha / 20 + h / 2
  }
  save(result, file = file.path(data_dir, paste0("heri_final_10_", h, "_n_5000_p_1000_", i, ".RData")))
}
run_reader <- function(label) {
  log <- file.path(work, paste0(label, ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(file.path(fixture, "run.R"), "simulations/annogk/figure2_result.R")),
                    stdout = log, stderr = log)
  list(status = status, log = paste(readLines(log, warn = FALSE), collapse = "\n"))
}
setwd(dirname(fixture))
ok <- run_reader("complete")
if (ok$status != 0) stop(ok$log)
out <- file.path(output_dir, "plots/main")
values <- read.delim(file.path(out, "AnnoGK_simu_strong_p1000_values.tsv"))
stopifnot(nrow(values) == 168L, !file.exists("Rplots.pdf"))
slots <- c(Knockoff = 1, AnnoKn = 3, GhostKnockoff = 4, AnnoGK = 6)
h <- as.numeric(sub(".*= ", "", values$heritability))
expected <- 50.5 / 1000 + unname(slots[values$method]) / 100 + values$alpha / 10 + h
expected[values$measure == "FDR"] <- expected[values$measure == "FDR"] / 2
stopifnot(max(abs(values$value - expected)) < 1e-12)
inputs <- read.delim(file.path(out, "AnnoGK_simu_strong_p1000_inputs.tsv"))
stopifnot(nrow(inputs) == 300L, all(nchar(inputs$md5) == 32L))
pdf <- file.path(out, "AnnoGK_simu_strong_p1000.pdf")
stopifnot(file.info(pdf)$size > 1000)
# Invalid inputs must fail before replacing any previously generated output.
old_hash <- tools::md5sum(pdf)
bad_file <- file.path(data_dir, "heri_final_10_0.05_n_5000_p_1000_1.RData")
load(bad_file)
result$power1[1] <- NA_real_
save(result, file = bad_file)
bad <- run_reader("invalid")
stopifnot(bad$status != 0, grepl("Invalid power1", bad$log, fixed = TRUE),
          identical(tools::md5sum(pdf), old_hash))
unlink(bad_file)
missing <- run_reader("missing")
stopifnot(missing$status != 0, grepl("Missing 1:", missing$log, fixed = TRUE),
          identical(tools::md5sum(pdf), old_hash))
cat("PASS: relocated standalone reader, 168 expected means, PDF/manifest export, no default plot device, invalid and missing input rejection.\n")
cat("Test files:", work, "\n")
