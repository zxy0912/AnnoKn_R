if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop(
    "Usage: Rscript map_genecards_height_top500_to_hg19.R [reference_dir]",
    call. = FALSE
  )
}
script_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", script_arguments, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[[1L]]))
} else normalizePath("map_genecards_height_top500_to_hg19.R")
analysis_dir <- repro_path("height_validation_workspace")
reference_dir <- if (length(args)) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else file.path(analysis_dir, "reference")

genecards_file <- file.path(
  reference_dir, "genecards_height_top500.csv"
)
gtf_file <- file.path(reference_dir, "gencode.v19.annotation.gtf.gz")
hgnc_file <- file.path(
  reference_dir, "hgnc_complete_set_2026-08-30.tsv"
)
hgnc_source_file <- file.path(
  reference_dir, "hgnc_complete_set_2026-08-30.source.txt"
)
required_files <- c(
  genecards_file, gtf_file, hgnc_file, hgnc_source_file
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required file(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", "),
    call. = FALSE
  )
}

.md5 <- function(path) unname(as.character(tools::md5sum(path)))
.collapse <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x)) paste(x, collapse = ";") else NA_character_
}
.clean_pipe_field <- function(value) {
  value <- gsub('"', "", ifelse(is.na(value), "", value), fixed = TRUE)
  if (!nzchar(value)) return(character())
  tokens <- trimws(strsplit(value, "|", fixed = TRUE)[[1L]])
  tokens[nzchar(tokens)]
}
.extract_attribute <- function(values, key) {
  pattern <- paste0("(?:^|; )", key, " \\\"([^\\\"]+)\\\"")
  matches <- regexec(pattern, values, perl = TRUE)
  pieces <- regmatches(values, matches)
  vapply(pieces, function(piece) {
    if (length(piece) >= 2L) piece[[2L]] else NA_character_
  }, character(1L))
}
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
.write_csv <- function(object, target_file) {
  temporary_file <- paste0(target_file, ".tmp.", Sys.getpid())
  write.csv(object, temporary_file, row.names = FALSE, na = "")
  .atomic_install(temporary_file, target_file)
}
.write_lines <- function(lines, target_file) {
  temporary_file <- paste0(target_file, ".tmp.", Sys.getpid())
  writeLines(lines, temporary_file)
  .atomic_install(temporary_file, target_file)
}
.write_bed4 <- function(object, target_file) {
  temporary_file <- paste0(target_file, ".tmp.", Sys.getpid())
  write.table(
    object, temporary_file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, na = ""
  )
  .atomic_install(temporary_file, target_file)
}

# Validate the frozen coordinate source before reading its features.
gtf_connection <- gzfile(gtf_file, open = "rt")
gtf_header <- readLines(gtf_connection, n = 4L, warn = FALSE)
close(gtf_connection)
gtf_header_text <- paste(gtf_header, collapse = " | ")
gtf_checks <- c(
  header_GRCh37 = grepl("GRCh37", gtf_header_text, fixed = TRUE),
  header_version_19 = grepl("version 19", gtf_header_text, fixed = TRUE),
  header_ensembl_74 = grepl("Ensembl 74", gtf_header_text, fixed = TRUE),
  official_release19_md5 = identical(
    .md5(gtf_file), "bd83e28270e595d3bde6bfcb21c9748f"
  )
)
if (!all(gtf_checks)) {
  stop(
    "GENCODE v19 build/checksum validation failed: ",
    paste(names(gtf_checks)[!gtf_checks], collapse = ", "),
    call. = FALSE
  )
}
if (!identical(
  .md5(hgnc_file), "acc6344dd46eb2942627a16e0d735fb8"
)) {
  stop("The frozen HGNC snapshot checksum changed.", call. = FALSE)
}

genecards <- read.csv(
  genecards_file, check.names = FALSE, stringsAsFactors = FALSE
)
expected_columns <- c(
  "Rank", "Symbol", "Name", "Type", "Relevance Score", "Knowledge"
)
input_checks <- c(
  columns = identical(names(genecards), expected_columns),
  rows = nrow(genecards) == 500L,
  ranks = identical(as.integer(genecards$Rank), 1:500),
  unique_symbols = !anyDuplicated(genecards$Symbol),
  symbols_present = all(nzchar(trimws(genecards$Symbol))),
  types = all(genecards$Type %in%
                c("Protein Coding", "RNA Gene", "Functional Element"))
)
if (!all(input_checks)) {
  stop(
    "GeneCards input validation failed: ",
    paste(names(input_checks)[!input_checks], collapse = ", "),
    call. = FALSE
  )
}

# Restrict the comprehensive GTF to gene features.  Coordinates in GTF are
# 1-based and closed; no transcript/exon union is substituted for the gene
# span supplied by GENCODE.
gtf_command <- sprintf(
  "gzip -cd %s | awk -F '\\t' '$3 == \"gene\" {print}'",
  shQuote(gtf_file)
)
gtf_raw <- read.delim(
  pipe(gtf_command), header = FALSE, sep = "\t", quote = "",
  comment.char = "", stringsAsFactors = FALSE
)
if (ncol(gtf_raw) != 9L || !nrow(gtf_raw)) {
  stop("Could not recover GENCODE gene features.", call. = FALSE)
}
gencode <- data.frame(
  seqname = gtf_raw$V1,
  source = gtf_raw$V2,
  feature = gtf_raw$V3,
  start_1based = as.integer(gtf_raw$V4),
  end_1based = as.integer(gtf_raw$V5),
  strand = gtf_raw$V7,
  gene_id_versioned = .extract_attribute(gtf_raw$V9, "gene_id"),
  gencode_v19_gene_name = .extract_attribute(gtf_raw$V9, "gene_name"),
  gencode_v19_gene_type = .extract_attribute(gtf_raw$V9, "gene_type"),
  gencode_v19_gene_status = .extract_attribute(
    gtf_raw$V9, "gene_status"
  ),
  stringsAsFactors = FALSE
)
gencode$ensembl_gene_id <- sub(
  "\\..*$", "", gencode$gene_id_versioned
)
gencode_checks <- c(
  features = all(gencode$feature == "gene"),
  finite_coordinates = all(is.finite(c(
    gencode$start_1based, gencode$end_1based
  ))),
  ordered_coordinates = all(
    gencode$start_1based <= gencode$end_1based
  ),
  identifiers = all(nzchar(gencode$gene_id_versioned)) &&
    all(nzchar(gencode$gencode_v19_gene_name)),
  reference_chromosomes = all(gencode$seqname %in%
    c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
)
if (!all(gencode_checks)) {
  stop(
    "GENCODE gene-feature validation failed: ",
    paste(names(gencode_checks)[!gencode_checks], collapse = ", "),
    call. = FALSE
  )
}
gencode_by_name <- split(
  seq_len(nrow(gencode)), gencode$gencode_v19_gene_name
)
gencode_by_ensembl <- split(
  seq_len(nrow(gencode)), gencode$ensembl_gene_id
)

hgnc <- read.delim(
  hgnc_file, check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE
)
hgnc_required <- c(
  "hgnc_id", "symbol", "name", "locus_type", "status",
  "location", "alias_symbol", "prev_symbol", "ensembl_gene_id"
)
if (any(!hgnc_required %in% names(hgnc))) {
  stop("The HGNC snapshot has an unexpected schema.", call. = FALSE)
}
hgnc_by_symbol <- split(seq_len(nrow(hgnc)), hgnc$symbol)
token_rows <- list()
token_index <- 0L
for (row_index in seq_len(nrow(hgnc))) {
  tokens <- unique(c(
    .clean_pipe_field(hgnc$prev_symbol[[row_index]]),
    .clean_pipe_field(hgnc$alias_symbol[[row_index]])
  ))
  if (length(tokens)) {
    for (token in tokens) {
      token_index <- token_index + 1L
      token_rows[[token_index]] <- data.frame(
        token = token, row = row_index, stringsAsFactors = FALSE
      )
    }
  }
}
hgnc_token_table <- do.call(rbind, token_rows)
hgnc_by_previous_or_alias <- split(
  hgnc_token_table$row, hgnc_token_table$token
)

.hgnc_rows <- function(symbol) {
  exact <- hgnc_by_symbol[[symbol]]
  if (length(exact)) return(list(rows = exact, type = "approved_symbol"))
  alias <- unique(hgnc_by_previous_or_alias[[symbol]])
  list(rows = alias, type = if (length(alias)) {
    "previous_or_alias_symbol"
  } else NA_character_)
}

detail_parts <- list()
detail_index <- 0L
for (input_index in seq_len(nrow(genecards))) {
  input <- genecards[input_index, , drop = FALSE]
  symbol <- input$Symbol[[1L]]
  exact_gene_rows <- gencode_by_name[[symbol]]
  hgnc_match <- .hgnc_rows(symbol)
  hgnc_rows <- hgnc_match$rows
  if (length(hgnc_rows) > 1L) {
    approved <- hgnc_rows[hgnc$status[hgnc_rows] == "Approved"]
    if (length(approved) == 1L) hgnc_rows <- approved
  }
  if (length(hgnc_rows) > 1L) {
    stop("Ambiguous HGNC mapping for symbol: ", symbol, call. = FALSE)
  }
  hgnc_row <- if (length(hgnc_rows)) hgnc_rows[[1L]] else NA_integer_

  gene_rows <- exact_gene_rows
  match_method <- if (length(gene_rows)) {
    "exact_gencode_v19_gene_name"
  } else NA_character_
  if (!length(gene_rows) && input$Type[[1L]] != "Functional Element" &&
      is.finite(hgnc_row)) {
    ensembl_ids <- .clean_pipe_field(
      hgnc$ensembl_gene_id[[hgnc_row]]
    )
    gene_rows <- unique(unlist(
      gencode_by_ensembl[ensembl_ids], use.names = FALSE
    ))
    if (length(gene_rows)) {
      match_method <- "hgnc_approved_symbol_to_ensembl_id"
    }
    if (!length(gene_rows)) {
      candidate_names <- unique(c(
        hgnc$symbol[[hgnc_row]],
        .clean_pipe_field(hgnc$prev_symbol[[hgnc_row]]),
        .clean_pipe_field(hgnc$alias_symbol[[hgnc_row]])
      ))
      gene_rows <- unique(unlist(
        gencode_by_name[candidate_names], use.names = FALSE
      ))
      if (length(gene_rows)) {
        match_method <- "hgnc_previous_or_alias_to_gencode_name"
      }
    }
  }

  if (!length(gene_rows)) {
    detail_index <- detail_index + 1L
    detail_parts[[detail_index]] <- data.frame(
      genecards_rank = as.integer(input$Rank),
      genecards_symbol = symbol,
      genecards_name = input$Name,
      genecards_type = input$Type,
      genecards_relevance_score = input[["Relevance Score"]],
      genecards_knowledge = input$Knowledge,
      mapping_status = if (input$Type[[1L]] == "Functional Element") {
        "not_a_gene_functional_element"
      } else "unmapped_gene_or_rna",
      match_method = NA_character_,
      hgnc_id = if (is.finite(hgnc_row)) hgnc$hgnc_id[hgnc_row] else NA,
      hgnc_approved_symbol = if (is.finite(hgnc_row)) {
        hgnc$symbol[hgnc_row]
      } else NA,
      hgnc_locus_type = if (is.finite(hgnc_row)) {
        hgnc$locus_type[hgnc_row]
      } else NA,
      hgnc_location = if (is.finite(hgnc_row)) {
        hgnc$location[hgnc_row]
      } else NA,
      gene_id_versioned = NA_character_,
      ensembl_gene_id = NA_character_,
      gencode_v19_gene_name = NA_character_,
      gencode_v19_gene_type = NA_character_,
      gencode_v19_gene_status = NA_character_,
      gencode_source = NA_character_,
      seqname = NA_character_, chromosome = NA_character_,
      start_1based = NA_integer_, end_1based = NA_integer_,
      width_bp = NA_integer_, strand = NA_character_,
      autosomal_1_22 = FALSE,
      coordinate_system = "GRCh37/hg19",
      interval_convention = "1-based closed",
      stringsAsFactors = FALSE
    )
    next
  }

  for (gene_row in gene_rows) {
    detail_index <- detail_index + 1L
    chromosome <- sub("^chr", "", gencode$seqname[[gene_row]])
    detail_parts[[detail_index]] <- data.frame(
      genecards_rank = as.integer(input$Rank),
      genecards_symbol = symbol,
      genecards_name = input$Name,
      genecards_type = input$Type,
      genecards_relevance_score = input[["Relevance Score"]],
      genecards_knowledge = input$Knowledge,
      mapping_status = if (identical(
        match_method, "exact_gencode_v19_gene_name"
      )) "mapped_exact_symbol" else "mapped_via_hgnc_ensembl",
      match_method = match_method,
      hgnc_id = if (is.finite(hgnc_row)) hgnc$hgnc_id[hgnc_row] else NA,
      hgnc_approved_symbol = if (is.finite(hgnc_row)) {
        hgnc$symbol[hgnc_row]
      } else NA,
      hgnc_locus_type = if (is.finite(hgnc_row)) {
        hgnc$locus_type[hgnc_row]
      } else NA,
      hgnc_location = if (is.finite(hgnc_row)) {
        hgnc$location[hgnc_row]
      } else NA,
      gene_id_versioned = gencode$gene_id_versioned[[gene_row]],
      ensembl_gene_id = gencode$ensembl_gene_id[[gene_row]],
      gencode_v19_gene_name =
        gencode$gencode_v19_gene_name[[gene_row]],
      gencode_v19_gene_type =
        gencode$gencode_v19_gene_type[[gene_row]],
      gencode_v19_gene_status =
        gencode$gencode_v19_gene_status[[gene_row]],
      gencode_source = gencode$source[[gene_row]],
      seqname = gencode$seqname[[gene_row]],
      chromosome = chromosome,
      start_1based = gencode$start_1based[[gene_row]],
      end_1based = gencode$end_1based[[gene_row]],
      width_bp = gencode$end_1based[[gene_row]] -
        gencode$start_1based[[gene_row]] + 1L,
      strand = gencode$strand[[gene_row]],
      autosomal_1_22 = chromosome %in% as.character(1:22),
      coordinate_system = "GRCh37/hg19",
      interval_convention = "1-based closed",
      stringsAsFactors = FALSE
    )
  }
}
detailed <- do.call(rbind, detail_parts)
rownames(detailed) <- NULL

# Collapse only overlapping GENCODE models for the same input symbol on the
# same chromosome.  Distinct chromosomes/loci (for example SHOX in PAR1) are
# retained as separate rows.  The detailed table always preserves every
# original GENCODE gene feature and Ensembl ID.
mapped <- detailed[grepl("^mapped_", detailed$mapping_status), , drop = FALSE]
collapsed_parts <- list()
collapsed_index <- 0L
mapped_groups <- split(
  seq_len(nrow(mapped)),
  interaction(mapped$genecards_rank, mapped$seqname, drop = TRUE)
)
for (group_indices in mapped_groups) {
  group <- mapped[group_indices, , drop = FALSE]
  group <- group[order(group$start_1based, group$end_1based), ]
  running_end <- -Inf
  locus_group <- integer(nrow(group))
  locus_number <- 0L
  for (row_position in seq_len(nrow(group))) {
    if (group$start_1based[[row_position]] > running_end) {
      locus_number <- locus_number + 1L
    }
    locus_group[[row_position]] <- locus_number
    running_end <- max(running_end, group$end_1based[[row_position]])
  }
  for (local_locus in unique(locus_group)) {
    part <- group[locus_group == local_locus, , drop = FALSE]
    collapsed_index <- collapsed_index + 1L
    collapsed_parts[[collapsed_index]] <- data.frame(
      genecards_rank = part$genecards_rank[[1L]],
      genecards_symbol = part$genecards_symbol[[1L]],
      genecards_name = part$genecards_name[[1L]],
      genecards_type = part$genecards_type[[1L]],
      genecards_relevance_score = part$genecards_relevance_score[[1L]],
      genecards_knowledge = part$genecards_knowledge[[1L]],
      mapping_status = .collapse(part$mapping_status),
      match_method = .collapse(part$match_method),
      hgnc_id = .collapse(part$hgnc_id),
      hgnc_approved_symbol = .collapse(part$hgnc_approved_symbol),
      hgnc_locus_type = .collapse(part$hgnc_locus_type),
      hgnc_location = .collapse(part$hgnc_location),
      gene_id_versioned = .collapse(part$gene_id_versioned),
      ensembl_gene_id = .collapse(part$ensembl_gene_id),
      gencode_v19_gene_name = .collapse(part$gencode_v19_gene_name),
      gencode_v19_gene_type = .collapse(part$gencode_v19_gene_type),
      seqname = part$seqname[[1L]],
      chromosome = part$chromosome[[1L]],
      start_1based = min(part$start_1based),
      end_1based = max(part$end_1based),
      width_bp = max(part$end_1based) - min(part$start_1based) + 1L,
      strand = .collapse(part$strand),
      gencode_model_count = nrow(part),
      autosomal_1_22 = part$autosomal_1_22[[1L]],
      coordinate_system = "GRCh37/hg19",
      interval_convention = "1-based closed",
      stringsAsFactors = FALSE
    )
  }
}
collapsed_mapped <- do.call(rbind, collapsed_parts)
unmapped <- detailed[!grepl("^mapped_", detailed$mapping_status), , drop = FALSE]
unmapped_collapsed <- data.frame(
  unmapped[, c(
    "genecards_rank", "genecards_symbol", "genecards_name",
    "genecards_type", "genecards_relevance_score",
    "genecards_knowledge", "mapping_status", "match_method",
    "hgnc_id", "hgnc_approved_symbol", "hgnc_locus_type",
    "hgnc_location"
  )],
  gene_id_versioned = NA_character_, ensembl_gene_id = NA_character_,
  gencode_v19_gene_name = NA_character_,
  gencode_v19_gene_type = NA_character_,
  seqname = NA_character_, chromosome = NA_character_,
  start_1based = NA_integer_, end_1based = NA_integer_,
  width_bp = NA_integer_, strand = NA_character_,
  gencode_model_count = 0L, autosomal_1_22 = FALSE,
  coordinate_system = "GRCh37/hg19",
  interval_convention = "1-based closed",
  stringsAsFactors = FALSE
)
positions <- rbind(collapsed_mapped, unmapped_collapsed)
locus_counts <- table(
  collapsed_mapped$genecards_rank
)
positions$locus_index <- NA_integer_
positions$locus_count <- 0L
for (rank_value in unique(collapsed_mapped$genecards_rank)) {
  indices <- which(positions$genecards_rank == rank_value &
                     !is.na(positions$start_1based))
  indices <- indices[order(
    match(positions$chromosome[indices], c(as.character(1:22), "X", "Y", "M")),
    positions$start_1based[indices]
  )]
  positions$locus_index[indices] <- seq_along(indices)
  positions$locus_count[indices] <- length(indices)
}
positions <- positions[order(
  positions$genecards_rank,
  ifelse(is.na(positions$locus_index), Inf, positions$locus_index)
), ]
rownames(positions) <- NULL
detailed <- detailed[order(
  detailed$genecards_rank,
  match(detailed$chromosome, c(as.character(1:22), "X", "Y", "M")),
  detailed$start_1based
), ]
rownames(detailed) <- NULL

nonfunctional <- genecards$Type != "Functional Element"
mapped_input_ranks <- unique(mapped$genecards_rank)
exact_input_ranks <- unique(mapped$genecards_rank[
  mapped$match_method == "exact_gencode_v19_gene_name"
])
hgnc_input_ranks <- unique(mapped$genecards_rank[
  mapped$match_method != "exact_gencode_v19_gene_name"
])
autosomal_input_ranks <- unique(
  collapsed_mapped$genecards_rank[collapsed_mapped$autosomal_1_22]
)
validation <- c(
  gtf_checks, input_checks, gencode_checks,
  all_gene_or_rna_inputs_mapped = setequal(
    genecards$Rank[nonfunctional], mapped_input_ranks
  ),
  no_functional_element_as_gene = !any(
    genecards$Rank[!nonfunctional] %in% mapped_input_ranks
  ),
  mapped_coordinates_valid = all(
    is.finite(c(mapped$start_1based, mapped$end_1based))
  ) && all(mapped$start_1based <= mapped$end_1based),
  output_contains_every_input = setequal(
    positions$genecards_rank, genecards$Rank
  ),
  output_rank_order = !is.unsorted(positions$genecards_rank),
  exact_plus_hgnc_partition = setequal(
    union(exact_input_ranks, hgnc_input_ranks), mapped_input_ranks
  ) && !length(intersect(exact_input_ranks, hgnc_input_ranks))
)
if (!all(validation)) {
  stop(
    "Final mapping validation failed: ",
    paste(names(validation)[!validation], collapse = ", "),
    call. = FALSE
  )
}

summary <- data.frame(
  metric = c(
    "genecards_input_rows", "unique_input_symbols",
    "protein_coding_inputs", "rna_gene_inputs",
    "functional_element_inputs", "gene_or_rna_inputs",
    "gene_or_rna_inputs_mapped", "exact_gencode_symbol_mappings",
    "hgnc_ensembl_rescue_mappings", "unmapped_gene_or_rna_inputs",
    "functional_elements_not_mapped_as_genes",
    "detailed_mapped_gencode_gene_models", "collapsed_mapped_gene_loci",
    "inputs_with_an_autosomal_locus", "nonautosomal_only_gene_inputs",
    "inputs_with_multiple_gene_loci", "final_position_rows",
    "validation_checks", "validation_checks_passed"
  ),
  value = as.character(c(
    nrow(genecards), length(unique(genecards$Symbol)),
    sum(genecards$Type == "Protein Coding"),
    sum(genecards$Type == "RNA Gene"),
    sum(genecards$Type == "Functional Element"), sum(nonfunctional),
    length(mapped_input_ranks), length(exact_input_ranks),
    length(hgnc_input_ranks),
    sum(nonfunctional) - length(mapped_input_ranks),
    nrow(unmapped), nrow(mapped), nrow(collapsed_mapped),
    length(autosomal_input_ranks),
    length(setdiff(mapped_input_ranks, autosomal_input_ranks)),
    sum(locus_counts > 1L), nrow(positions),
    length(validation), sum(validation)
  )),
  stringsAsFactors = FALSE
)

positions_file <- file.path(
  reference_dir, "genecards_height_top500_hg19_gene_positions.csv"
)
detailed_file <- file.path(
  reference_dir,
  "genecards_height_top500_hg19_gene_models_detailed.csv"
)
unmapped_file <- file.path(
  reference_dir,
  "genecards_height_top500_hg19_non_gene_or_unmapped.csv"
)
summary_file <- file.path(
  reference_dir, "genecards_height_top500_hg19_mapping_summary.csv"
)
metadata_file <- file.path(
  reference_dir, "genecards_height_top500_hg19_mapping_metadata.txt"
)
bed_file <- file.path(
  reference_dir, "genecards_height_top500_hg19_gene_positions.bed"
)
autosomal_bed_file <- file.path(
  reference_dir,
  "genecards_height_top500_hg19_autosomal_gene_positions.bed"
)
completion_file <- file.path(
  reference_dir, "GENECARDS_HEIGHT_TOP500_HG19_MAPPING_COMPLETE.txt"
)

bed <- data.frame(
  chrom = collapsed_mapped$seqname,
  start_0based = collapsed_mapped$start_1based - 1L,
  end_0based_exclusive = collapsed_mapped$end_1based,
  name = paste0(
    "rank", collapsed_mapped$genecards_rank, "|",
    collapsed_mapped$genecards_symbol, "|locus",
    ave(
      collapsed_mapped$genecards_rank,
      collapsed_mapped$genecards_rank,
      FUN = seq_along
    )
  ),
  stringsAsFactors = FALSE
)
bed <- bed[order(
  match(sub("^chr", "", bed$chrom), c(as.character(1:22), "X", "Y", "M")),
  bed$start_0based, bed$end_0based_exclusive
), ]
autosomal_bed <- bed[sub("^chr", "", bed$chrom) %in%
                        as.character(1:22), , drop = FALSE]

.write_csv(positions, positions_file)
.write_csv(detailed, detailed_file)
.write_csv(unmapped, unmapped_file)
.write_csv(summary, summary_file)
.write_bed4(bed, bed_file)
.write_bed4(autosomal_bed, autosomal_bed_file)

metadata <- c(
  "experiment=GeneCards height top-500 symbol to GRCh37/hg19 gene-span mapping",
  "mapping_snapshot_date=2026-08-30",
  paste0("implementation=", script_file),
  paste0("implementation_md5=", .md5(script_file)),
  paste0("genecards_file=", genecards_file),
  paste0("genecards_md5=", .md5(genecards_file)),
  paste0("gencode_file=", gtf_file),
  paste0("gencode_md5=", .md5(gtf_file)),
  paste0("gencode_header=", gtf_header_text),
  "gencode_official_release=https://www.gencodegenes.org/human/release_19.html",
  "gencode_official_md5_source=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/MD5SUMS",
  paste0("hgnc_file=", hgnc_file),
  paste0("hgnc_md5=", .md5(hgnc_file)),
  paste0("hgnc_source_metadata=", hgnc_source_file),
  "coordinate_source=GENCODE v19 comprehensive gene features",
  "coordinate_system=GRCh37.p13/hg19",
  "csv_interval_convention=1-based closed",
  "bed_interval_convention=0-based half-open",
  paste(
    "mapping_rule=exact GENCODE v19 gene_name first; otherwise current",
    "HGNC approved/previous/alias nomenclature to stable Ensembl gene ID;",
    "functional elements are not mapped as genes"
  ),
  paste(
    "collapse_rule=merge only overlapping GENCODE gene models for the same",
    "GeneCards rank and chromosome; retain distinct loci/chromosomes"
  ),
  "height_analysis_chromosomes=1-22 only",
  paste(
    "build_audit_warning=height.R directly intersects GIANT Position with",
    "PLINK bim positions and performs no liftOver; the inaccessible bychr_*",
    "reference panel must still be independently verified as GRCh37"
  ),
  paste0("validation=", paste(names(validation), validation, sep = ":", collapse = ";"))
)
.write_lines(metadata, metadata_file)

output_files <- c(
  positions_file, detailed_file, unmapped_file, summary_file,
  metadata_file, bed_file, autosomal_bed_file
)
completion <- c(
  "status=complete",
  "coordinate_system=GRCh37.p13/hg19",
  "csv_interval_convention=1-based closed",
  "bed_interval_convention=0-based half-open",
  paste0("input_rows=", nrow(genecards)),
  paste0("gene_or_rna_inputs_mapped=", length(mapped_input_ranks)),
  paste0("functional_elements_not_mapped_as_genes=", nrow(unmapped)),
  paste0("collapsed_mapped_gene_loci=", nrow(collapsed_mapped)),
  paste0("validation_checks=", length(validation)),
  "all_validations_passed=true",
  paste0("implementation_md5=", .md5(script_file)),
  paste0("output_md5=", paste(
    basename(output_files), vapply(output_files, .md5, character(1L)),
    sep = ":", collapse = ";"
  ))
)
.write_lines(completion, completion_file)
cat("Completed GeneCards top-500 GRCh37/hg19 gene-position mapping.\n")
