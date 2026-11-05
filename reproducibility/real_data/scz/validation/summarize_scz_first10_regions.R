if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
# Fixed first ten seeds, chosen before inspecting intersections.
input_dir <- repro_path("scz_validation_workspace", "result")
output_dir <- file.path(dirname(input_dir), "first10_seed_region_stability_q0p1")
seeds <- seq(100L, 1000L, 100L)
methods <- c("GhostKnockoff", "AnnoGK", "Transfer", "AdaKn")
records <- list()
for (ancestry in c("EUR", "EAS")) {
  other <- if (ancestry == "EUR") "EAS_AFR_LAT" else "EUR_AFR_LAT"
  for (seed in seeds) for (chr in setdiff(1:22, 21L)) {
    path <- file.path(input_dir, sprintf(
      "%s_result_adakn_10_median_r2_%s_chr_%d_M_1_%d.RData", ancestry, other, chr, seed))
    e <- new.env()
    load(path, envir = e)
    r <- e$result
    coords <- t(vapply(r$risk_region, function(x) as.numeric(unlist(x)), numeric(2)))
    hits <- lapply(r[c("q_gk", "q_gkanno", "q_transfer")], function(qs)
      vapply(qs, function(q) {stopifnot(!anyNA(q)); any(q <= 0.1)}, logical(1)))
    hits[[4L]] <- lengths(r$rej_adakn) > 0L
    names(hits) <- methods
    records[[length(records) + 1L]] <- data.frame(ancestry, seed, chr,
      start = coords[, 1], end = coords[, 2],
      region_id = paste(chr, coords[, 1], coords[, 2], sep = "_"), hits)
  }
}
stopifnot(length(records) == 420L)
x <- do.call(rbind, records)
stopifnot(!anyDuplicated(x[c("ancestry", "seed", "region_id")]))
groups <- c("ancestry", "chr", "start", "end", "region_id")
observed <- aggregate(list(n_seeds = x$seed), x[groups], length)
stopifnot(all(observed$n_seeds == 10L))
freq <- aggregate(x[methods], x[groups], sum)
summary <- do.call(rbind, lapply(c("EUR", "EAS"), function(a) {
  r <- freq[freq$ancestry == a, ]
  do.call(rbind, lapply(methods, function(m) data.frame(ancestry = a, method = m,
    union_regions = sum(r[[m]] >= 1L), all10_regions = sum(r[[m]] == 10L),
    maximum_seed_recurrence = max(r[[m]]))))
}))
# Cross-check raw-file reconstruction against the previous 100-seed extraction.
previous <- read.delim(file.path(dirname(input_dir),
  "hundred_seed_region_stability_q0p1", "region_by_seed_detections.tsv"))
previous <- previous[previous$seed %in% seeds, ]
keys <- function(z) paste(z$ancestry, z$seed, z$region_id, sep = ":")
idx <- match(keys(x), keys(previous))
stopifnot(nrow(x) == nrow(previous), !anyNA(idx))
for (m in methods) stopifnot(identical(x[[m]], previous[[m]][idx]))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary, file.path(output_dir, "method_region_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(freq, file.path(output_dir, "region_seed_recurrence.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(list(seeds = seeds, M = 1L, q_threshold = 0.1,
  chromosomes = setdiff(1:22, 21L), n_files = length(records),
  rule = "All ten seeds detect at least one SNP in the same coordinate-defined region",
  AdaKn_rule = "Nonempty stored rejection-index vector",
  raw_files_match_previous_extraction = TRUE, session = sessionInfo()),
  file = file.path(output_dir, "metadata.txt"))
print(summary, row.names = FALSE)
