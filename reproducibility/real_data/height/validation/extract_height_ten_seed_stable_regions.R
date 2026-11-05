if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1L) args[[1L]] else repro_path("height_validation_workspace", "annot_pvalus/result")
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  repro_path("height_validation_workspace", "height_ten_seed_region_stability_q0p1")
}
q_target <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.1

ancestry <- "AFR"
seeds <- seq.int(100L, 1000L, by = 100L)
chromosomes <- seq_len(22L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result_path <- function(chromosome, seed) {
  file.path(
    result_dir,
    sprintf(
      "%s_result_final_10_median_r2_EUR_chr_%d_M_1_%d.RData",
      ancestry,
      chromosome,
      seed
    )
  )
}

expected_paths <- unlist(lapply(chromosomes, function(chromosome) {
  vapply(seeds, function(seed) result_path(chromosome, seed), character(1L))
}))
missing_paths <- expected_paths[!file.exists(expected_paths)]
if (length(missing_paths)) {
  stop("Missing expected result files: ", paste(missing_paths, collapse = ", "))
}
if (length(expected_paths) != 220L) stop("Expected exactly 220 AFR result files")

detected <- function(q_values) {
  q_values <- as.numeric(q_values)
  any(is.finite(q_values) & q_values <= q_target)
}

minimum_q <- function(q_values) {
  q_values <- as.numeric(q_values)
  q_values <- q_values[is.finite(q_values)]
  if (length(q_values)) min(q_values) else NA_real_
}

rows <- vector("list", length(expected_paths))
file_index <- 0L

for (chromosome in chromosomes) {
  for (seed in seeds) {
    file_index <- file_index + 1L
    path <- result_path(chromosome, seed)
    loaded <- new.env(parent = emptyenv())
    load(path, envir = loaded)
    if (!exists("result", envir = loaded, inherits = FALSE)) {
      stop("Object 'result' is missing from ", path)
    }
    result <- loaded$result
    if (as.integer(result$chr) != chromosome) {
      stop("Chromosome mismatch in ", path)
    }
    n_regions <- length(result$risk_region)
    if (
      length(result$q_gk) != n_regions ||
        length(result$q_gkanno) != n_regions
    ) {
      stop("Region/q-value list length mismatch in ", path)
    }

    region_rows <- vector("list", n_regions)
    for (region_index in seq_len(n_regions)) {
      interval <- as.numeric(unlist(result$risk_region[[region_index]]))
      if (length(interval) != 2L || any(!is.finite(interval))) {
        stop("Invalid risk region in ", path, " at index ", region_index)
      }
      start <- interval[[1L]]
      end <- interval[[2L]]
      if (start > end) stop("Reversed risk region in ", path)

      region_rows[[region_index]] <- data.frame(
        ancestry = ancestry,
        chromosome = chromosome,
        start = start,
        end = end,
        region_id = sprintf("chr%d:%s-%s", chromosome, start, end),
        seed = seed,
        ghost_detected = detected(result$q_gk[[region_index]]),
        annogk_detected = detected(result$q_gkanno[[region_index]]),
        ghost_min_q = minimum_q(result$q_gk[[region_index]]),
        annogk_min_q = minimum_q(result$q_gkanno[[region_index]]),
        source_file = basename(path),
        stringsAsFactors = FALSE
      )
    }
    file_rows <- if (n_regions) do.call(rbind, region_rows) else data.frame()
    if (nrow(file_rows) && anyDuplicated(file_rows$region_id)) {
      stop("Duplicated region within one seed/file: ", path)
    }
    rows[[file_index]] <- file_rows
  }
}

detection_by_seed <- do.call(rbind, rows)
detection_by_seed <- detection_by_seed[
  order(
    detection_by_seed$chromosome,
    detection_by_seed$start,
    detection_by_seed$end,
    detection_by_seed$seed
  ),
]

region_keys <- unique(detection_by_seed[, c(
  "ancestry", "chromosome", "start", "end", "region_id"
)])
region_keys <- region_keys[
  order(region_keys$chromosome, region_keys$start, region_keys$end),
]

summarize_region <- function(region_id) {
  x <- detection_by_seed[detection_by_seed$region_id == region_id, ]
  data.frame(
    region_id = region_id,
    seeds_with_region_object = length(unique(x$seed)),
    ghost_n_seeds = length(unique(x$seed[x$ghost_detected])),
    annogk_n_seeds = length(unique(x$seed[x$annogk_detected])),
    ghost_worst_detected_min_q = if (any(x$ghost_detected)) {
      max(x$ghost_min_q[x$ghost_detected], na.rm = TRUE)
    } else {
      NA_real_
    },
    annogk_worst_detected_min_q = if (any(x$annogk_detected)) {
      max(x$annogk_min_q[x$annogk_detected], na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

region_summaries <- do.call(rbind, lapply(region_keys$region_id, summarize_region))
recurrence <- merge(region_keys, region_summaries, by = "region_id", sort = FALSE)
recurrence <- recurrence[
  order(recurrence$chromosome, recurrence$start, recurrence$end),
]
recurrence$ghost_all_10 <- recurrence$ghost_n_seeds == length(seeds)
recurrence$annogk_all_10 <- recurrence$annogk_n_seeds == length(seeds)

thresholds <- 2L:10L
recurrence_summary <- data.frame(
  ancestry = ancestry,
  at_least_k_seeds = thresholds,
  ghost_regions = vapply(
    thresholds,
    function(k) sum(recurrence$ghost_n_seeds >= k),
    integer(1L)
  ),
  annogk_regions = vapply(
    thresholds,
    function(k) sum(recurrence$annogk_n_seeds >= k),
    integer(1L)
  ),
  stringsAsFactors = FALSE
)

expected_ghost <- c(245L, 183L, 154L, 119L, 97L, 81L, 56L, 36L, 15L)
expected_annogk <- c(260L, 207L, 173L, 143L, 116L, 86L, 64L, 41L, 21L)
if (
  !identical(recurrence_summary$ghost_regions, expected_ghost) ||
    !identical(recurrence_summary$annogk_regions, expected_annogk)
) {
  print(recurrence_summary, row.names = FALSE)
  stop("Recomputed recurrence counts do not match the frozen handoff table")
}

all_10 <- recurrence[recurrence$ghost_all_10 | recurrence$annogk_all_10, ]
all_10$all_10_class <- ifelse(
  all_10$ghost_all_10 & all_10$annogk_all_10,
  "Both",
  ifelse(all_10$ghost_all_10, "GhostKnockoff_only", "AnnoGK_only")
)
all_10 <- all_10[
  order(all_10$chromosome, all_10$start, all_10$end),
]
ghost_all_10 <- recurrence[recurrence$ghost_all_10, ]
annogk_all_10 <- recurrence[recurrence$annogk_all_10, ]

write.table(
  detection_by_seed,
  file.path(output_dir, "AFR_region_detection_by_seed.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  recurrence,
  file.path(output_dir, "AFR_region_recurrence.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  recurrence_summary,
  file.path(output_dir, "AFR_recurrence_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  all_10,
  file.path(output_dir, "AFR_all10_union_regions.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  ghost_all_10,
  file.path(output_dir, "AFR_GhostKnockoff_all10_regions.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  annogk_all_10,
  file.path(output_dir, "AFR_AnnoGK_all10_regions.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Validated result files:", length(expected_paths), "\n")
cat("GhostKnockoff regions detected in all 10 seeds:", nrow(ghost_all_10), "\n")
cat("AnnoGK regions detected in all 10 seeds:", nrow(annogk_all_10), "\n")
cat("Union:", nrow(all_10), "\n")
print(all_10[, c(
  "region_id", "ghost_n_seeds", "annogk_n_seeds", "all_10_class"
)], row.names = FALSE)
