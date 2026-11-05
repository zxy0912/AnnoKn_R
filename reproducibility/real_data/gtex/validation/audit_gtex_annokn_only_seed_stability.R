if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
setwd(repro_path("gtex_validation_workspace"))
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
old_result_dir <- if (length(args) >= 1L) args[[1L]] else "result"
new_result_dir <- if (length(args) >= 2L) args[[2L]] else "new_seed/result"
candidate_file <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  paste0(
    "gtex_annokn_only_eGenes_q01_initial_seed_PET/",
    "annokn_only_egenes_with_independent_eqtl_support.tsv"
  )
}
output_dir <- if (length(args) >= 4L) {
  args[[4L]]
} else {
  "gtex_annokn_only_eGenes_q01_initial_seed_PET"
}
q_target <- if (length(args) >= 5L) as.numeric(args[[5L]]) else 0.1

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stable_gene_id <- function(x) sub("[.].*$", "", as.character(x))

knockoff_plus_threshold <- function(w, q) {
  candidates <- sort(unique(abs(w[is.finite(w) & w != 0])))
  if (!length(candidates)) return(Inf)
  ratios <- vapply(
    candidates,
    function(t) (1 + sum(w <= -t, na.rm = TRUE)) /
      max(1, sum(w >= t, na.rm = TRUE)),
    numeric(1)
  )
  passing <- which(ratios <= q)
  if (!length(passing)) Inf else candidates[min(passing)]
}

selected_indices <- function(record, method, q) {
  if (method == "knockoff") {
    w <- as.numeric(record$ghostknockoff)
  } else if (method == "annokn") {
    beta <- as.numeric(record$knockoff_anno$beta)
    if (length(beta) %% 2L != 0L) stop("Odd AnnoKn beta length")
    p <- length(beta) %/% 2L
    w <- abs(beta[seq_len(p)]) - abs(beta[p + seq_len(p)])
  } else {
    stop("Unknown method: ", method)
  }
  tau <- knockoff_plus_threshold(w, q)
  list(index = which(w >= tau), tau = tau, max_w = max(w))
}

selection_category <- function(is_knockoff, is_annokn) {
  if (is_annokn && !is_knockoff) return("AnnoKn_only")
  if (!is_annokn && is_knockoff) return("Knockoff_only")
  if (is_annokn && is_knockoff) return("Both")
  "Neither"
}

load_all_result <- function(path) {
  loaded <- new.env(parent = emptyenv())
  load(path, envir = loaded)
  if (!exists("all_result", envir = loaded, inherits = FALSE)) {
    stop("all_result is missing from ", path)
  }
  loaded$all_result
}

candidates <- read.delim(candidate_file, check.names = FALSE)
if (nrow(candidates) != 612L) {
  warning("Expected 612 old-seed AnnoKn-only pairs; found ", nrow(candidates))
}

required_columns <- c(
  "tissue", "source_file", "ensembl_gene_id", "hgnc_symbol",
  "selection_category", "n_selected_snps_knockoff",
  "n_selected_snps_annokn"
)
missing_columns <- setdiff(required_columns, names(candidates))
if (length(missing_columns)) {
  stop("Candidate file is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (any(candidates$selection_category != "AnnoKn_only")) {
  stop("Candidate input contains rows that are not old-seed AnnoKn-only")
}

audit_rows <- vector("list", nrow(candidates))
output_index <- 0L

for (filename in unique(candidates$source_file)) {
  old_path <- file.path(old_result_dir, filename)
  new_path <- file.path(new_result_dir, filename)
  if (!file.exists(old_path)) stop("Missing old result: ", old_path)
  if (!file.exists(new_path)) stop("Missing new result: ", new_path)

  old_results <- load_all_result(old_path)
  new_results <- load_all_result(new_path)
  old_ids <- vapply(
    old_results,
    function(record) stable_gene_id(record$gene_info$ensembl_gene_id[[1L]]),
    character(1)
  )
  new_ids <- vapply(
    new_results,
    function(record) stable_gene_id(record$gene_info$ensembl_gene_id[[1L]]),
    character(1)
  )
  if (!identical(old_ids, new_ids)) {
    stop("Old/new gene order mismatch in ", filename)
  }

  file_candidates <- candidates[candidates$source_file == filename, , drop = FALSE]
  for (candidate_index in seq_len(nrow(file_candidates))) {
    candidate <- file_candidates[candidate_index, , drop = FALSE]
    gene_id <- stable_gene_id(candidate$ensembl_gene_id)
    record_index <- match(gene_id, old_ids)
    if (is.na(record_index)) {
      stop("Candidate ", gene_id, " is absent from ", filename)
    }

    old_record <- old_results[[record_index]]
    new_record <- new_results[[record_index]]
    old_knockoff <- selected_indices(old_record, "knockoff", q_target)
    old_annokn <- selected_indices(old_record, "annokn", q_target)
    new_knockoff <- selected_indices(new_record, "knockoff", q_target)
    new_annokn <- selected_indices(new_record, "annokn", q_target)

    old_category <- selection_category(
      length(old_knockoff$index) > 0L,
      length(old_annokn$index) > 0L
    )
    if (
      old_category != "AnnoKn_only" ||
      length(old_knockoff$index) != candidate$n_selected_snps_knockoff ||
      length(old_annokn$index) != candidate$n_selected_snps_annokn
    ) {
      stop("Old-seed validation failed for ", candidate$tissue, " / ", gene_id)
    }

    new_category <- selection_category(
      length(new_knockoff$index) > 0L,
      length(new_annokn$index) > 0L
    )
    old_annokn_indices <- old_annokn$index
    new_annokn_indices <- new_annokn$index
    shared_indices <- intersect(old_annokn_indices, new_annokn_indices)
    union_indices <- union(old_annokn_indices, new_annokn_indices)

    output_index <- output_index + 1L
    audit_rows[[output_index]] <- cbind(
      candidate,
      data.frame(
        new_seed_category = new_category,
        new_seed_annokn_egene = length(new_annokn_indices) > 0L,
        new_seed_knockoff_egene = length(new_knockoff$index) > 0L,
        new_seed_n_selected_snps_annokn = length(new_annokn_indices),
        new_seed_n_selected_snps_knockoff = length(new_knockoff$index),
        new_seed_tau_annokn = new_annokn$tau,
        new_seed_tau_knockoff = new_knockoff$tau,
        old_new_annokn_shared_snp_columns = length(shared_indices),
        old_new_annokn_snp_jaccard = if (length(union_indices)) {
          length(shared_indices) / length(union_indices)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    )
  }
}

audit <- do.call(rbind, audit_rows)
audit <- audit[order(audit$tissue, audit$chromosome, audit$ensembl_gene_id), ]

tissues <- sort(unique(audit$tissue))
summary_rows <- lapply(tissues, function(tissue) {
  x <- audit[audit$tissue == tissue, , drop = FALSE]
  category_counts <- table(factor(
    x$new_seed_category,
    levels = c("AnnoKn_only", "Both", "Knockoff_only", "Neither")
  ))
  retained <- sum(x$new_seed_annokn_egene)
  data.frame(
    tissue = tissue,
    old_seed_annokn_only = nrow(x),
    retained_as_new_seed_annokn_egene = retained,
    annokn_egene_retention_fraction = retained / nrow(x),
    remains_annokn_only = unname(category_counts[["AnnoKn_only"]]),
    becomes_both = unname(category_counts[["Both"]]),
    becomes_knockoff_only = unname(category_counts[["Knockoff_only"]]),
    becomes_neither = unname(category_counts[["Neither"]]),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
overall_counts <- table(factor(
  audit$new_seed_category,
  levels = c("AnnoKn_only", "Both", "Knockoff_only", "Neither")
))
overall <- data.frame(
  tissue = "Overall",
  old_seed_annokn_only = nrow(audit),
  retained_as_new_seed_annokn_egene = sum(audit$new_seed_annokn_egene),
  annokn_egene_retention_fraction = mean(audit$new_seed_annokn_egene),
  remains_annokn_only = unname(overall_counts[["AnnoKn_only"]]),
  becomes_both = unname(overall_counts[["Both"]]),
  becomes_knockoff_only = unname(overall_counts[["Knockoff_only"]]),
  becomes_neither = unname(overall_counts[["Neither"]]),
  stringsAsFactors = FALSE
)
summary_table <- rbind(summary_table, overall)

write.table(
  audit,
  file.path(output_dir, "old_annokn_only_new_seed_status.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  summary_table,
  file.path(output_dir, "old_annokn_only_new_seed_summary_by_tissue.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Validated old-seed candidates:", nrow(audit), "\n")
cat(
  "Retained as AnnoKn eGenes under new seed:",
  sum(audit$new_seed_annokn_egene),
  sprintf("(%.1f%%)\n", 100 * mean(audit$new_seed_annokn_egene))
)
print(summary_table, row.names = FALSE)
