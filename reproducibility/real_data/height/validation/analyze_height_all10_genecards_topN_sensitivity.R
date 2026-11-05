if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

# Sensitivity analysis for the GeneCards Height rank cutoff.  For each
# Top-N set (N = 100, 200, ..., 1000), summarize the maximum GeneCards
# relevance score per stable AFR Height region and then compare the region
# classes and the two method-level sets.  Regions without a Top-N hit are
# retained in hit-count denominators but are not assigned a score of zero.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
analysis_dir <- repro_path("height_validation_workspace")
stability_dir <- file.path(
  analysis_dir, "height_ten_seed_region_stability_q0p1"
)

region_file <- file.path(stability_dir, "AFR_all10_union_regions.tsv")
long_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_overlap_long.tsv"
)
overlap_marker <- file.path(
  stability_dir,
  "AFR_ALL10_GENECARDS_HEIGHT_TOP1000_OVERLAP_COMPLETE.txt"
)

class_output_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_topN_class_sensitivity.tsv"
)
method_output_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_topN_method_sensitivity.tsv"
)
comparison_output_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_topN_mean_comparison.tsv"
)
strict_comparison_output_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_topN_mean_comparison_strict.tsv"
)
metadata_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_topN_sensitivity_metadata.txt"
)
completion_file <- file.path(
  stability_dir,
  "AFR_ALL10_GENECARDS_HEIGHT_TOPN_SENSITIVITY_COMPLETE.txt"
)

required_files <- c(region_file, long_file, overlap_marker)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", "),
    call. = FALSE
  )
}

.md5 <- function(path) unname(as.character(tools::md5sum(path)))

.atomic_install <- function(temporary_file, target_file) {
  if (file.exists(target_file)) {
    if (identical(.md5(temporary_file), .md5(target_file))) {
      unlink(temporary_file)
      return(invisible(target_file))
    }
    stop(
      "Refusing to replace a nonidentical existing output: ", target_file,
      call. = FALSE
    )
  }
  if (!file.rename(temporary_file, target_file)) {
    unlink(temporary_file)
    stop("Could not atomically install: ", target_file, call. = FALSE)
  }
  invisible(target_file)
}

.write_tsv <- function(object, target_file) {
  temporary_file <- paste0(target_file, ".tmp.", Sys.getpid())
  write.table(
    object, temporary_file, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = ""
  )
  .atomic_install(temporary_file, target_file)
}

.write_lines <- function(lines, target_file) {
  temporary_file <- paste0(target_file, ".tmp.", Sys.getpid())
  writeLines(lines, temporary_file)
  .atomic_install(temporary_file, target_file)
}

.conditional_summary <- function(region_subset) {
  scores <- region_subset$region_max_relevance_score
  scores <- scores[is.finite(scores)]
  data.frame(
    total_regions = nrow(region_subset),
    regions_with_topN_hit = length(scores),
    regions_without_topN_hit = nrow(region_subset) - length(scores),
    hit_fraction = length(scores) / nrow(region_subset),
    conditional_mean_region_max_score = if (length(scores)) {
      mean(scores)
    } else {
      NA_real_
    },
    conditional_median_region_max_score = if (length(scores)) {
      median(scores)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

regions <- read.delim(region_file, check.names = FALSE)
long <- read.delim(long_file, check.names = FALSE)

required_region_columns <- c(
  "region_id", "all_10_class", "ghost_all_10", "annogk_all_10"
)
required_long_columns <- c(
  "region_id", "all_10_class", "genecards_rank", "relevance_score",
  "fully_contained_in_region"
)
if (!all(required_region_columns %in% names(regions))) {
  stop("Unexpected stable-region table schema.", call. = FALSE)
}
if (!all(required_long_columns %in% names(long))) {
  stop("Unexpected top-1000 overlap-long table schema.", call. = FALSE)
}

long$genecards_rank <- as.integer(long$genecards_rank)
long$relevance_score <- as.numeric(long$relevance_score)
if (!is.logical(long$fully_contained_in_region)) {
  long$fully_contained_in_region <- as.logical(long$fully_contained_in_region)
}

expected_class_counts <- c(
  Both = 9L, GhostKnockoff_only = 6L, AnnoGK_only = 12L
)
actual_class_counts <- table(factor(
  regions$all_10_class, levels = names(expected_class_counts)
))
input_checks <- c(
  regions_27 = nrow(regions) == 27L,
  unique_regions = !anyDuplicated(regions$region_id),
  expected_class_counts = identical(
    unname(as.integer(actual_class_counts)), unname(expected_class_counts)
  ),
  long_regions_known = all(long$region_id %in% regions$region_id),
  long_region_gene_pairs_unique = !anyDuplicated(paste(
    long$region_id, long$genecards_rank, sep = "|"
  )),
  valid_ranks = all(long$genecards_rank %in% 1:1000),
  finite_scores = all(is.finite(long$relevance_score)),
  valid_containment_flags = !anyNA(long$fully_contained_in_region),
  class_labels_agree = all(
    long$all_10_class == regions$all_10_class[
      match(long$region_id, regions$region_id)
    ]
  ),
  upstream_complete = any(
    readLines(overlap_marker, warn = FALSE) == "all_validations_passed=true"
  )
)
if (!all(input_checks)) {
  stop(
    "Input validation failed: ",
    paste(names(input_checks)[!input_checks], collapse = ", "),
    call. = FALSE
  )
}

cutoffs <- seq.int(100L, 1000L, by = 100L)
overlap_rules <- c("any_gene_body_overlap", "fully_contained_gene")
class_levels <- c("Both", "GhostKnockoff_only", "AnnoGK_only")
class_labels <- c(
  Both = "Both", GhostKnockoff_only = "Ghost-only",
  AnnoGK_only = "AnnoGK-only"
)
method_levels <- c("GhostKnockoff", "AnnoGK")

.region_scores <- function(top_n, overlap_rule) {
  hits <- long[long$genecards_rank <= top_n, , drop = FALSE]
  if (identical(overlap_rule, "fully_contained_gene")) {
    hits <- hits[hits$fully_contained_in_region, , drop = FALSE]
  }

  maxima <- if (nrow(hits)) {
    aggregate(relevance_score ~ region_id, hits, max)
  } else {
    data.frame(region_id = character(), relevance_score = numeric())
  }
  names(maxima)[names(maxima) == "relevance_score"] <-
    "region_max_relevance_score"

  result <- merge(
    regions[, required_region_columns, drop = FALSE], maxima,
    by = "region_id", all.x = TRUE, sort = FALSE
  )
  result <- result[match(regions$region_id, result$region_id), ]
  rownames(result) <- NULL
  result
}

class_parts <- list()
method_parts <- list()
class_index <- 0L
method_index <- 0L

for (overlap_rule in overlap_rules) {
  for (top_n in cutoffs) {
    scores <- .region_scores(top_n, overlap_rule)

    for (class_name in class_levels) {
      class_index <- class_index + 1L
      one_summary <- .conditional_summary(
        scores[scores$all_10_class == class_name, , drop = FALSE]
      )
      class_parts[[class_index]] <- data.frame(
        top_n = top_n,
        overlap_rule = overlap_rule,
        all_10_class = class_name,
        display_class = unname(class_labels[class_name]),
        one_summary,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    for (method_name in method_levels) {
      method_index <- method_index + 1L
      selected <- if (identical(method_name, "GhostKnockoff")) {
        scores$ghost_all_10
      } else {
        scores$annogk_all_10
      }
      one_summary <- .conditional_summary(scores[selected, , drop = FALSE])
      method_parts[[method_index]] <- data.frame(
        top_n = top_n,
        overlap_rule = overlap_rule,
        method = method_name,
        included_classes = if (identical(method_name, "GhostKnockoff")) {
          "Both + Ghost-only"
        } else {
          "Both + AnnoGK-only"
        },
        one_summary,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
}

class_summary <- do.call(rbind, class_parts)
method_summary <- do.call(rbind, method_parts)
rownames(class_summary) <- NULL
rownames(method_summary) <- NULL

.wide_comparison <- function(overlap_rule) {
  class_subset <- class_summary[
    class_summary$overlap_rule == overlap_rule, , drop = FALSE
  ]
  method_subset <- method_summary[
    method_summary$overlap_rule == overlap_rule, , drop = FALSE
  ]
  output <- data.frame(top_n = cutoffs)

  for (class_name in class_levels) {
    rows <- class_subset[
      class_subset$all_10_class == class_name,
      , drop = FALSE
    ]
    rows <- rows[match(cutoffs, rows$top_n), ]
    prefix <- c(
      Both = "both", GhostKnockoff_only = "ghost_only",
      AnnoGK_only = "annogk_only"
    )[[class_name]]
    output[[paste0(prefix, "_hit_regions")]] <- rows$regions_with_topN_hit
    output[[paste0(prefix, "_conditional_mean")]] <-
      rows$conditional_mean_region_max_score
  }

  for (method_name in method_levels) {
    rows <- method_subset[
      method_subset$method == method_name,
      , drop = FALSE
    ]
    rows <- rows[match(cutoffs, rows$top_n), ]
    prefix <- if (identical(method_name, "GhostKnockoff")) {
      "ghost_including_both"
    } else {
      "annogk_including_both"
    }
    output[[paste0(prefix, "_hit_regions")]] <- rows$regions_with_topN_hit
    output[[paste0(prefix, "_conditional_mean")]] <-
      rows$conditional_mean_region_max_score
  }

  output$ghost_only_minus_annogk_only_mean <-
    output$ghost_only_conditional_mean - output$annogk_only_conditional_mean
  output$ghost_minus_annogk_including_both_mean <-
    output$ghost_including_both_conditional_mean -
    output$annogk_including_both_conditional_mean
  output
}

comparison <- .wide_comparison("any_gene_body_overlap")
strict_comparison <- .wide_comparison("fully_contained_gene")

.close <- function(observed, expected, tolerance = 1e-8) {
  isTRUE(all.equal(as.numeric(observed), as.numeric(expected), tolerance = tolerance))
}

row_200 <- comparison[comparison$top_n == 200L, ]
row_500 <- comparison[comparison$top_n == 500L, ]
row_1000 <- comparison[comparison$top_n == 1000L, ]
strict_1000 <- strict_comparison[strict_comparison$top_n == 1000L, ]

output_checks <- c(
  class_rows_60 = nrow(class_summary) == 60L,
  method_rows_40 = nrow(method_summary) == 40L,
  comparison_rows_10 = nrow(comparison) == 10L,
  strict_comparison_rows_10 = nrow(strict_comparison) == 10L,
  hit_counts_nondecreasing = all(vapply(
    split(
      class_summary$regions_with_topN_hit,
      interaction(class_summary$overlap_rule, class_summary$all_10_class)
    ),
    function(values) all(diff(values) >= 0L), logical(1)
  )) && all(vapply(
    split(
      method_summary$regions_with_topN_hit,
      interaction(method_summary$overlap_rule, method_summary$method)
    ),
    function(values) all(diff(values) >= 0L), logical(1)
  )),
  top200_both_mean = .close(row_200$both_conditional_mean, 92.95),
  top200_ghost_only_mean = .close(
    row_200$ghost_only_conditional_mean, 51.7666666666667
  ),
  top200_annogk_only_mean = .close(
    row_200$annogk_only_conditional_mean, 60.2
  ),
  top500_both_mean = .close(row_500$both_conditional_mean, 75),
  top500_ghost_only_mean = .close(
    row_500$ghost_only_conditional_mean, 49.75
  ),
  top500_annogk_only_mean = .close(
    row_500$annogk_only_conditional_mean, 45.3
  ),
  top1000_both_mean = .close(
    row_1000$both_conditional_mean, 51.1666666666667
  ),
  top1000_ghost_only_mean = .close(
    row_1000$ghost_only_conditional_mean, 41.6
  ),
  top1000_annogk_only_mean = .close(
    row_1000$annogk_only_conditional_mean, 40.7875
  ),
  strict_top1000_ghost_only_mean = .close(
    strict_1000$ghost_only_conditional_mean, 39.7166666666667
  ),
  strict_top1000_annogk_only_mean = .close(
    strict_1000$annogk_only_conditional_mean, 40.7875
  )
)
if (!all(output_checks)) {
  stop(
    "Output validation failed: ",
    paste(names(output_checks)[!output_checks], collapse = ", "),
    call. = FALSE
  )
}

.write_tsv(class_summary, class_output_file)
.write_tsv(method_summary, method_output_file)
.write_tsv(comparison, comparison_output_file)
.write_tsv(strict_comparison, strict_comparison_output_file)

metadata_lines <- c(
  "Analysis: GeneCards Height Top-N sensitivity for AFR all-10-seed regions",
  "Date: 2026-08-30",
  paste0("Region input: ", normalizePath(region_file)),
  paste0("Region input MD5: ", .md5(region_file)),
  paste0("Top-1000 overlap input: ", normalizePath(long_file)),
  paste0("Top-1000 overlap input MD5: ", .md5(long_file)),
  "Cutoffs: 100,200,300,400,500,600,700,800,900,1000",
  "Region score: maximum GeneCards relevance score among Top-N genes in the region",
  paste0(
    "Main overlap rule: gene_start <= region_end AND gene_end >= region_start ",
    "(GRCh37/hg19, 1-based closed)"
  ),
  "Strict sensitivity rule: gene_start >= region_start AND gene_end <= region_end",
  paste0(
    "Missing-data rule: regions without a Top-N hit remain missing and are ",
    "excluded from conditional score means; they are never scored as zero"
  ),
  "GhostKnockoff method set: Both + GhostKnockoff_only",
  "AnnoGK method set: Both + AnnoGK_only",
  paste0("Input validations passed: ", sum(input_checks), "/", length(input_checks)),
  paste0("Output validations passed: ", sum(output_checks), "/", length(output_checks)),
  "all_validations_passed=true"
)
.write_lines(metadata_lines, metadata_file)

completion_lines <- c(
  "GENECARDS HEIGHT TOP-N SENSITIVITY ANALYSIS COMPLETE",
  "all_validations_passed=true",
  paste0("input_validations_passed=", sum(input_checks), "/", length(input_checks)),
  paste0("output_validations_passed=", sum(output_checks), "/", length(output_checks)),
  paste0("class_output_md5=", .md5(class_output_file)),
  paste0("method_output_md5=", .md5(method_output_file)),
  paste0("comparison_output_md5=", .md5(comparison_output_file)),
  paste0("strict_comparison_output_md5=", .md5(strict_comparison_output_file)),
  paste0("metadata_md5=", .md5(metadata_file)),
  paste0("implementation_md5=", .md5(script_file))
)
.write_lines(completion_lines, completion_file)

cat("GeneCards Height Top-N sensitivity analysis complete.\n")
cat("Input validations:", sum(input_checks), "/", length(input_checks), "passed.\n")
cat("Output validations:", sum(output_checks), "/", length(output_checks), "passed.\n")
cat("Main comparison:\n")
print(comparison, row.names = FALSE)
cat("Strict-containment comparison:\n")
print(strict_comparison, row.names = FALSE)
