if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
setwd(repro_path("gtex_validation_workspace"))
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1L) args[[1L]] else "result"
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "gtex_annokn_only_eGenes_q01_initial_seed_PET"
}
q_target <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.1
validation_mode <- if (length(args) >= 4L) tolower(args[[4L]]) else "full"
if (!validation_mode %in% c("full", "structure")) {
  stop("Validation mode must be either 'full' or 'structure'.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

collapse_values <- function(x) {
  if (!length(x)) "" else paste(x, collapse = ";")
}

files <- sort(list.files(result_dir, pattern = "[.]RData$", full.names = TRUE))
if (!length(files)) stop("No RData files found in ", result_dir)

gene_rows <- vector("list", 0L)
validation_rows <- vector("list", length(files))
row_index <- 0L

for (file_index in seq_along(files)) {
  file <- files[[file_index]]
  filename <- basename(file)
  match_parts <- regexec(
    "^chr_([0-9]+)_ct_(.+)_type_([0-9]+)[.]RData$",
    filename
  )
  parsed <- regmatches(filename, match_parts)[[1L]]
  if (length(parsed) != 4L) stop("Unexpected filename: ", filename)

  chromosome <- as.integer(parsed[[2L]])
  tissue_key <- parsed[[3L]]
  tissue <- gsub("_", " ", tissue_key, fixed = TRUE)

  loaded <- new.env(parent = emptyenv())
  load(file, envir = loaded)
  if (!exists("all_result", envir = loaded, inherits = FALSE)) {
    stop("all_result is missing from ", file)
  }
  all_result <- loaded$all_result

  n_knockoff <- 0L
  n_annokn <- 0L

  for (gene_index in seq_along(all_result)) {
    record <- all_result[[gene_index]]
    info <- record$gene_info[1L, , drop = FALSE]

    w_knockoff <- as.numeric(record$ghostknockoff)
    tau_knockoff <- knockoff_plus_threshold(w_knockoff, q_target)
    selected_knockoff <- which(w_knockoff >= tau_knockoff)

    beta <- as.numeric(record$knockoff_anno$beta)
    if (length(beta) %% 2L != 0L) {
      stop("Odd AnnoKn beta length in ", filename, ", gene ", gene_index)
    }
    p <- length(beta) %/% 2L
    w_annokn <- abs(beta[seq_len(p)]) - abs(beta[p + seq_len(p)])
    tau_annokn <- knockoff_plus_threshold(w_annokn, q_target)
    selected_annokn <- which(w_annokn >= tau_annokn)

    is_knockoff <- length(selected_knockoff) > 0L
    is_annokn <- length(selected_annokn) > 0L
    n_knockoff <- n_knockoff + is_knockoff
    n_annokn <- n_annokn + is_annokn

    category <- if (is_annokn && !is_knockoff) {
      "AnnoKn_only"
    } else if (!is_annokn && is_knockoff) {
      "Knockoff_only"
    } else if (is_annokn && is_knockoff) {
      "Both"
    } else {
      "Neither"
    }

    if (category == "Neither") next

    snp_names <- as.character(record$snp_name)
    annotation <- as.data.frame(record$annot)
    annotation_columns <- intersect(
      c("Promoter_UCSC", "Enhancer_Hoffman", "TFBS_ENCODE"),
      names(annotation)
    )

    annotation_hits <- setNames(rep(0L, length(annotation_columns)),
                                annotation_columns)
    if (length(selected_annokn) && length(annotation_columns)) {
      annotation_hits <- vapply(
        annotation[selected_annokn, annotation_columns, drop = FALSE],
        function(x) sum(as.numeric(x) != 0, na.rm = TRUE),
        integer(1)
      )
    }

    lambda <- rep(0, length(annotation_columns))
    names(lambda) <- annotation_columns
    raw_lambda <- unlist(record$knockoff_anno$lambda_s, use.names = FALSE)
    if (length(annotation_columns)) {
      annotation_matrix <- as.matrix(record$R)
      annotation_variances <- apply(
        annotation_matrix, 2L, var, na.rm = TRUE
      )
      variable_columns <- is.finite(annotation_variances) &
        annotation_variances != 0
      if (sum(variable_columns) == length(raw_lambda)) {
        lambda[variable_columns] <- raw_lambda
      } else if (length(raw_lambda) == length(lambda)) {
        lambda[] <- raw_lambda
      }
    }

    row_index <- row_index + 1L
    gene_rows[[row_index]] <- data.frame(
      tissue = tissue,
      tissue_key = tissue_key,
      chromosome = chromosome,
      ensembl_gene_id = as.character(info$ensembl_gene_id),
      hgnc_symbol = as.character(info$hgnc_symbol),
      start_position = as.integer(info$start_position),
      end_position = as.integer(info$end_position),
      selection_category = category,
      n_selected_snps_knockoff = length(selected_knockoff),
      n_selected_snps_annokn = length(selected_annokn),
      selected_snps_knockoff = collapse_values(snp_names[selected_knockoff]),
      selected_snps_annokn = collapse_values(snp_names[selected_annokn]),
      max_w_knockoff = if (length(w_knockoff)) max(w_knockoff) else NA_real_,
      max_w_annokn = if (length(w_annokn)) max(w_annokn) else NA_real_,
      tau_knockoff = tau_knockoff,
      tau_annokn = tau_annokn,
      selected_annokn_promoter_hits = unname(annotation_hits["Promoter_UCSC"]),
      selected_annokn_enhancer_hits = unname(annotation_hits["Enhancer_Hoffman"]),
      selected_annokn_tfbs_hits = unname(annotation_hits["TFBS_ENCODE"]),
      lambda_promoter = unname(lambda["Promoter_UCSC"]),
      lambda_enhancer = unname(lambda["Enhancer_Hoffman"]),
      lambda_tfbs = unname(lambda["TFBS_ENCODE"]),
      source_file = filename,
      stringsAsFactors = FALSE
    )
  }

  validation_rows[[file_index]] <- data.frame(
    file = filename,
    chr = chromosome,
    celltype = tissue_key,
    q = q_target,
    n_gene = length(all_result),
    n_eGene_knockoff = n_knockoff,
    n_eGene_AnnoKn = n_annokn,
    stringsAsFactors = FALSE
  )
}

selected <- do.call(rbind, gene_rows)
validation <- do.call(rbind, validation_rows)

annokn_only <- selected[selected$selection_category == "AnnoKn_only", ]
knockoff_only <- selected[selected$selection_category == "Knockoff_only", ]

count_levels <- c("AnnoKn_only", "Knockoff_only", "Both")
counts <- as.data.frame.matrix(table(
  factor(selected$tissue, levels = sort(unique(selected$tissue))),
  factor(selected$selection_category, levels = count_levels)
))
counts$tissue <- rownames(counts)
rownames(counts) <- NULL
counts <- counts[, c("tissue", count_levels)]
counts$n_eGene_knockoff <- counts$Knockoff_only + counts$Both
counts$n_eGene_AnnoKn <- counts$AnnoKn_only + counts$Both
counts$net_gain <- counts$n_eGene_AnnoKn - counts$n_eGene_knockoff

detail_reference <- read.csv(
  "sensitivity_annotation_settings_detail.csv",
  check.names = FALSE
)
expected <- detail_reference[
  detail_reference$condition == "full" & detail_reference$q == q_target,
  c("file", "n_gene", "n_eGene_knockoff", "n_eGene_AnnoKn")
]
validation_check <- merge(
  validation,
  expected,
  by = "file",
  suffixes = c("_recomputed", "_saved"),
  all = TRUE
)
validation_check$file_and_gene_structure_matches <- with(
  validation_check,
  !is.na(n_gene_recomputed) & !is.na(n_gene_saved) &
    n_gene_recomputed == n_gene_saved
)
validation_check$discovery_counts_match_reference <- with(
  validation_check,
  file_and_gene_structure_matches &
    n_eGene_knockoff_recomputed == n_eGene_knockoff_saved &
    n_eGene_AnnoKn_recomputed == n_eGene_AnnoKn_saved
)
validation_check$matches <- if (validation_mode == "full") {
  validation_check$discovery_counts_match_reference
} else {
  validation_check$file_and_gene_structure_matches
}

if (anyNA(validation_check$matches) || !all(validation_check$matches)) {
  write.table(
    validation_check,
    file.path(output_dir, "validation_by_file.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  stop(
    if (validation_mode == "full") {
      "Recomputed selection counts did not match saved reference counts."
    } else {
      "Recomputed file/gene structure did not match the reference analysis."
    }
  )
}

write.table(
  annokn_only,
  file.path(output_dir, "annokn_only_egenes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  knockoff_only,
  file.path(output_dir, "knockoff_only_egenes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  selected,
  gzfile(file.path(output_dir, "all_selected_egenes.tsv.gz")),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  counts,
  file.path(output_dir, "selection_counts_by_tissue.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  validation_check,
  file.path(output_dir, "validation_by_file.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Validated files:", nrow(validation_check), "\n")
cat("Validation mode:", validation_mode, "\n")
cat("AnnoKn-only eGenes:", nrow(annokn_only), "\n")
cat("Knockoff-only eGenes:", nrow(knockoff_only), "\n")
cat("Output directory:", normalizePath(output_dir), "\n")
print(counts, row.names = FALSE)
