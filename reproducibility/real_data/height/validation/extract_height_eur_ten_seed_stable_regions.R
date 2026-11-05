if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

# EUR counterpart of the frozen AFR extraction. Only M=1 and the ten
# seeds 100, 200, ..., 1000 are included. No model fitting is performed.
options(stringsAsFactors = FALSE)
script_file <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
))
analysis_dir <- repro_path("height_validation_workspace")
result_dir <- file.path(analysis_dir, "annot_pvalus", "result")
output_dir <- file.path(analysis_dir, "height_ten_seed_region_stability_q0p1")
seeds <- seq.int(100L, 1000L, by = 100L)
q_target <- 0.1

inputs <- expand.grid(seed = seeds, chromosome = seq_len(22L))
inputs$source_file <- sprintf(
  "EUR_result_final_10_median_r2_AFR_chr_%d_M_1_%d.RData",
  inputs$chromosome, inputs$seed
)
inputs$path <- file.path(result_dir, inputs$source_file)
stopifnot(nrow(inputs) == 220L, !anyDuplicated(inputs$path))
if (!all(file.exists(inputs$path))) {
  stop("Missing expected EUR inputs: ", paste(inputs$path[!file.exists(inputs$path)], collapse = ", "))
}

minimum_q <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

parts <- vector("list", nrow(inputs))
for (i in seq_len(nrow(inputs))) {
  loaded <- new.env(parent = emptyenv())
  load(inputs$path[i], envir = loaded)
  if (!exists("result", envir = loaded, inherits = FALSE)) stop("Missing result: ", inputs$path[i])
  result <- loaded$result
  stopifnot(length(result$chr) == 1L, as.integer(result$chr) == inputs$chromosome[i])
  n_regions <- length(result$risk_region)
  stopifnot(length(result$q_gk) == n_regions, length(result$q_gkanno) == n_regions)
  if (!n_regions) next
  intervals <- lapply(result$risk_region, function(x) {
    interval <- as.numeric(unlist(x))
    stopifnot(length(interval) == 2L, all(is.finite(interval)), interval[1] <= interval[2])
    interval
  })
  starts <- vapply(intervals, `[`, numeric(1L), 1L)
  ends <- vapply(intervals, `[`, numeric(1L), 2L)
  ghost_min_q <- vapply(result$q_gk, minimum_q, numeric(1L))
  annogk_min_q <- vapply(result$q_gkanno, minimum_q, numeric(1L))
  one <- data.frame(
    ancestry = "EUR", chromosome = inputs$chromosome[i], start = starts, end = ends,
    region_id = sprintf("chr%d:%.0f-%.0f", inputs$chromosome[i], starts, ends),
    seed = inputs$seed[i],
    ghost_detected = is.finite(ghost_min_q) & ghost_min_q <= q_target,
    annogk_detected = is.finite(annogk_min_q) & annogk_min_q <= q_target,
    ghost_min_q = ghost_min_q, annogk_min_q = annogk_min_q,
    source_file = inputs$source_file[i]
  )
  stopifnot(!anyDuplicated(one$region_id))
  parts[[i]] <- one
}
detection <- do.call(rbind, parts)
detection <- detection[order(detection$chromosome, detection$start, detection$end, detection$seed), ]
stopifnot(!anyDuplicated(paste(detection$region_id, detection$seed, sep = "|")))

summaries <- lapply(split(detection, detection$region_id), function(x) {
  worst <- function(q, selected) if (any(selected)) max(q[selected]) else NA_real_
  cbind(
    x[1L, c("region_id", "ancestry", "chromosome", "start", "end")],
    data.frame(
      seeds_with_region_object = length(unique(x$seed)),
      ghost_n_seeds = length(unique(x$seed[x$ghost_detected])),
      annogk_n_seeds = length(unique(x$seed[x$annogk_detected])),
      ghost_worst_detected_min_q = worst(x$ghost_min_q, x$ghost_detected),
      annogk_worst_detected_min_q = worst(x$annogk_min_q, x$annogk_detected)
    )
  )
})
recurrence <- do.call(rbind, summaries)
recurrence <- recurrence[order(recurrence$chromosome, recurrence$start, recurrence$end), ]
rownames(recurrence) <- NULL
recurrence$ghost_all_10 <- recurrence$ghost_n_seeds == length(seeds)
recurrence$annogk_all_10 <- recurrence$annogk_n_seeds == length(seeds)
curve <- data.frame(
  ancestry = "EUR", at_least_k_seeds = 2L:10L,
  ghost_regions = vapply(2L:10L, function(k) sum(recurrence$ghost_n_seeds >= k), integer(1L)),
  annogk_regions = vapply(2L:10L, function(k) sum(recurrence$annogk_n_seeds >= k), integer(1L))
)
# Validate against the independently saved August 20 handoff, not only k=10.
stopifnot(
  identical(curve$ghost_regions, c(459L, 375L, 325L, 281L, 245L, 211L, 187L, 158L, 99L)),
  identical(curve$annogk_regions, c(593L, 508L, 446L, 395L, 355L, 322L, 284L, 246L, 192L))
)

all10 <- recurrence[recurrence$ghost_all_10 | recurrence$annogk_all_10, ]
all10$all_10_class <- ifelse(
  all10$ghost_all_10 & all10$annogk_all_10, "Both",
  ifelse(all10$annogk_all_10, "AnnoGK_only", "GhostKnockoff_only")
)
inputs$md5 <- unname(tools::md5sum(inputs$path))
inputs$bytes <- file.info(inputs$path)$size
inputs$path <- NULL
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_frozen <- function(x, name) {
  target <- file.path(output_dir, name)
  temporary <- tempfile(pattern = paste0(name, "."), tmpdir = output_dir)
  write.table(x, temporary, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  if (file.exists(target)) {
    same <- unname(tools::md5sum(temporary)) == unname(tools::md5sum(target))
    unlink(temporary)
    if (!same) stop("Existing output differs: ", target)
  } else if (!file.rename(temporary, target)) {
    stop("Could not install: ", target)
  }
}
write_frozen(detection, "EUR_region_detection_by_seed.tsv")
write_frozen(recurrence, "EUR_region_recurrence.tsv")
write_frozen(curve, "EUR_recurrence_summary.tsv")
write_frozen(all10, "EUR_all10_union_regions.tsv")
write_frozen(all10[all10$annogk_all_10, ], "EUR_AnnoGK_all10_regions.tsv")
write_frozen(all10[all10$ghost_all_10, ], "EUR_GhostKnockoff_all10_regions.tsv")
write_frozen(inputs, "EUR_ten_seed_input_manifest.tsv")
metadata <- data.frame(
  key = c("ancestry", "M", "seeds", "q_threshold", "expected_input_files", "recurrence_handoff_matches",
          "region_identity", "detection_rule", "script_md5", "R_version"),
  value = c("EUR", "1", paste(seeds, collapse = ","), "0.1", "220", "TRUE",
            "chromosome,start,end", "at least one finite SNP q-value <= 0.1 within each seed",
            unname(tools::md5sum(script_file)), R.version.string)
)
write_frozen(metadata, "EUR_ten_seed_extraction_metadata.tsv")

cat("Validated EUR input files:", nrow(inputs), "\n")
cat("All recurrence thresholds reproduce the frozen handoff.\n")
print(curve, row.names = FALSE)
cat("All-ten region classes:\n")
print(table(all10$all_10_class))
cat("Union regions:", nrow(all10), "\n")
