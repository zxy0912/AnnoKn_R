if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

# Intersect the 27 AFR-target Height regions detected in all ten M=1 seeds by
# GhostKnockoff and/or AnnoGK with the frozen GeneCards Height top-1000 gene
# spans on GRCh37/hg19.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
analysis_dir <- repro_path("height_validation_workspace")
stability_dir <- file.path(
  analysis_dir, "height_ten_seed_region_stability_q0p1"
)
reference_dir <- file.path(analysis_dir, "reference")

region_file <- file.path(stability_dir, "AFR_all10_union_regions.tsv")
gene_file <- file.path(
  reference_dir, "genecards_height_top1000_hg19_gene_positions.csv"
)
gene_mapping_marker <- file.path(
  reference_dir, "GENECARDS_HEIGHT_TOP1000_HG19_MAPPING_COMPLETE.txt"
)
top500_long_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top500_overlap_long.tsv"
)
top500_overlap_marker <- file.path(
  stability_dir,
  "AFR_ALL10_GENECARDS_HEIGHT_TOP500_OVERLAP_COMPLETE.txt"
)

summary_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_overlap.tsv"
)
long_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_overlap_long.tsv"
)
class_summary_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_overlap_class_summary.tsv"
)
ranking_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_ranked.tsv"
)
metadata_file <- file.path(
  stability_dir,
  "AFR_all10_union_regions_genecards_height_top1000_overlap_metadata.txt"
)
completion_file <- file.path(
  stability_dir,
  "AFR_ALL10_GENECARDS_HEIGHT_TOP1000_OVERLAP_COMPLETE.txt"
)

required_files <- c(
  region_file, gene_file, gene_mapping_marker,
  top500_long_file, top500_overlap_marker
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", "),
    call. = FALSE
  )
}

.md5 <- function(path) unname(as.character(tools::md5sum(path)))

.normalize_table_for_regression <- function(object) {
  normalized <- as.data.frame(
    lapply(object, function(values) {
      values <- as.character(values)
      values[is.na(values)] <- ""
      values
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(normalized) <- NULL
  normalized
}

.collapse <- function(values, separator = ";") {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) paste(values, collapse = separator) else ""
}

.format_score <- function(values) {
  format(as.numeric(values), scientific = FALSE, trim = TRUE, digits = 15L)
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

regions <- read.delim(region_file, check.names = FALSE)
required_region_columns <- c(
  "region_id", "ancestry", "chromosome", "start", "end",
  "seeds_with_region_object", "ghost_n_seeds", "annogk_n_seeds",
  "ghost_all_10", "annogk_all_10", "all_10_class"
)
if (!all(required_region_columns %in% names(regions))) {
  stop("Unexpected all-10 region table schema.", call. = FALSE)
}

regions$chromosome <- as.integer(regions$chromosome)
regions$start <- as.integer(regions$start)
regions$end <- as.integer(regions$end)
expected_class_counts <- c(Both = 9L, AnnoGK_only = 12L, GhostKnockoff_only = 6L)
actual_class_counts <- table(factor(
  regions$all_10_class, levels = names(expected_class_counts)
))
reconstructed_region_id <- paste0(
  "chr", regions$chromosome, ":", regions$start, "-", regions$end
)

region_checks <- c(
  rows_27 = nrow(regions) == 27L,
  unique_region_ids = !anyDuplicated(regions$region_id),
  ancestry_AFR = all(regions$ancestry == "AFR"),
  chromosomes_1_22 = all(regions$chromosome %in% 1:22),
  valid_intervals = all(
    is.finite(regions$start) & is.finite(regions$end) &
      regions$start >= 1L & regions$start <= regions$end
  ),
  reconstructed_ids = identical(regions$region_id, reconstructed_region_id),
  region_objects_in_all_seeds = all(regions$seeds_with_region_object == 10L),
  expected_class_counts = identical(
    unname(as.integer(actual_class_counts)), unname(expected_class_counts)
  ),
  class_logic = all(
    (regions$all_10_class == "Both") ==
      (regions$ghost_all_10 & regions$annogk_all_10)
  ) && all(
    (regions$all_10_class == "GhostKnockoff_only") ==
      (regions$ghost_all_10 & !regions$annogk_all_10)
  ) && all(
    (regions$all_10_class == "AnnoGK_only") ==
      (!regions$ghost_all_10 & regions$annogk_all_10)
  )
)
if (!all(region_checks)) {
  stop(
    "Region validation failed: ",
    paste(names(region_checks)[!region_checks], collapse = ", "),
    call. = FALSE
  )
}

gene_positions <- read.csv(gene_file, check.names = FALSE)
required_gene_columns <- c(
  "genecards_rank", "genecards_symbol", "genecards_name", "genecards_type",
  "genecards_relevance_score", "genecards_knowledge", "mapping_status",
  "match_method", "hgnc_approved_symbol", "ensembl_gene_id",
  "gencode_v19_gene_name", "gencode_v19_gene_type", "seqname",
  "chromosome", "start_1based", "end_1based", "width_bp", "strand",
  "autosomal_1_22", "coordinate_system", "interval_convention",
  "locus_index", "locus_count"
)
if (!all(required_gene_columns %in% names(gene_positions))) {
  stop("Unexpected GeneCards gene-position table schema.", call. = FALSE)
}

mapped <- gene_positions[grepl("^mapped_", gene_positions$mapping_status), ]
mapped$start_1based <- as.integer(mapped$start_1based)
mapped$end_1based <- as.integer(mapped$end_1based)
mapped$genecards_rank <- as.integer(mapped$genecards_rank)
mapped$genecards_relevance_score <- as.numeric(
  mapped$genecards_relevance_score
)

gene_checks <- c(
  all_1000_inputs_represented = identical(
    sort(unique(gene_positions$genecards_rank)), 1:1000
  ),
  mapped_gene_rna_inputs_990 = length(unique(mapped$genecards_rank)) == 990L,
  mapped_loci_991 = nrow(mapped) == 991L,
  no_functional_elements = !any(mapped$genecards_type == "Functional Element"),
  hg19_coordinates = all(mapped$coordinate_system == "GRCh37/hg19"),
  closed_intervals = all(mapped$interval_convention == "1-based closed"),
  valid_gene_intervals = all(
    mapped$start_1based >= 1L &
      mapped$start_1based <= mapped$end_1based &
      mapped$width_bp == mapped$end_1based - mapped$start_1based + 1L
  ),
  finite_relevance_scores = all(is.finite(mapped$genecards_relevance_score)),
  mapping_completion_marker = any(
    readLines(gene_mapping_marker, warn = FALSE) ==
      "all_validations_passed=true"
  )
)
if (!all(gene_checks)) {
  stop(
    "Gene-position validation failed: ",
    paste(names(gene_checks)[!gene_checks], collapse = ", "),
    call. = FALSE
  )
}

# The all-10 Height regions are restricted to chromosomes 1--22. Filter the
# mapped GeneCards loci before coercing chromosome labels so chrX/chrY do not
# become NA and accidentally propagate through R's interval predicates.
mapped <- mapped[mapped$autosomal_1_22, , drop = FALSE]
mapped$chromosome <- as.integer(mapped$chromosome)
if (nrow(mapped) != 936L || anyNA(mapped$chromosome)) {
  stop("Unexpected autosomal GeneCards locus set.", call. = FALSE)
}

class_levels <- c("Both", "AnnoGK_only", "GhostKnockoff_only")
display_class <- c(
  Both = "Both", AnnoGK_only = "AnnoGK-only",
  GhostKnockoff_only = "Ghost-only"
)
regions$class_order <- match(regions$all_10_class, class_levels)
regions <- regions[order(regions$class_order, regions$chromosome, regions$start), ]
rownames(regions) <- NULL

region_columns <- c(
  "region_id", "ancestry", "chromosome", "start", "end",
  "all_10_class", "ghost_n_seeds", "annogk_n_seeds"
)
long_parts <- list()
long_index <- 0L
summary_parts <- vector("list", nrow(regions))

for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index, ]
  overlaps <-
    mapped$chromosome == region$chromosome &
    mapped$start_1based <= region$end &
    mapped$end_1based >= region$start
  genes <- mapped[overlaps, , drop = FALSE]
  genes <- genes[order(
    -genes$genecards_relevance_score,
    genes$genecards_rank,
    genes$genecards_symbol,
    genes$start_1based
  ), ]
  fully_contained <-
    genes$start_1based >= region$start & genes$end_1based <= region$end

  if (nrow(genes)) {
    long_index <- long_index + 1L
    overlap_start <- pmax(region$start, genes$start_1based)
    overlap_end <- pmin(region$end, genes$end_1based)
    overlap_bp <- overlap_end - overlap_start + 1L
    gene_span_bp <- genes$end_1based - genes$start_1based + 1L
    long_parts[[long_index]] <- data.frame(
      region[, region_columns, drop = FALSE],
      display_class = unname(display_class[region$all_10_class]),
      genecards_rank = genes$genecards_rank,
      genecards_symbol = genes$genecards_symbol,
      genecards_name = genes$genecards_name,
      genecards_type = genes$genecards_type,
      relevance_score = genes$genecards_relevance_score,
      knowledge = genes$genecards_knowledge,
      match_method = genes$match_method,
      hgnc_approved_symbol = genes$hgnc_approved_symbol,
      gencode_v19_gene_name = genes$gencode_v19_gene_name,
      ensembl_gene_id = genes$ensembl_gene_id,
      gene_start_1based = genes$start_1based,
      gene_end_1based = genes$end_1based,
      gene_strand = genes$strand,
      overlap_start_1based = overlap_start,
      overlap_end_1based = overlap_end,
      overlap_bp = overlap_bp,
      gene_span_bp = gene_span_bp,
      fraction_gene_span_overlapped = overlap_bp / gene_span_bp,
      fully_contained_in_region = fully_contained,
      coordinate_system = "GRCh37/hg19",
      interval_convention = "1-based closed",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = NULL
    )
  }

  contains <- nrow(genes) > 0L
  maximum_score <- if (contains) max(genes$genecards_relevance_score) else NA_real_
  maximum_rows <- if (contains) {
    genes$genecards_relevance_score == maximum_score
  } else logical()
  fully_contained_maximum_score <- if (any(fully_contained)) {
    max(genes$genecards_relevance_score[fully_contained])
  } else NA_real_
  fully_contained_maximum_rows <- if (any(fully_contained)) {
    fully_contained &
      genes$genecards_relevance_score == fully_contained_maximum_score
  } else rep(FALSE, nrow(genes))

  summary_parts[[region_index]] <- data.frame(
    region[, region_columns, drop = FALSE],
    display_class = unname(display_class[region$all_10_class]),
    contains_genecards_top1000_gene = contains,
    top1000_gene_overlap_status = if (contains) {
      "contains_top1000_gene"
    } else {
      "does_not_contain_top1000_gene"
    },
    n_genecards_top1000_genes = nrow(genes),
    n_fully_contained_top1000_genes = sum(fully_contained),
    boundary_overlap_only_top1000_gene_symbols = .collapse(
      genes$genecards_symbol[!fully_contained]
    ),
    top1000_gene_symbols_by_relevance = .collapse(genes$genecards_symbol),
    top1000_gene_ranks_by_relevance = .collapse(genes$genecards_rank),
    top1000_relevance_scores_by_relevance = .collapse(
      if (contains) .format_score(genes$genecards_relevance_score) else character()
    ),
    top1000_gene_rank_score = .collapse(
      if (contains) {
        paste0(
          genes$genecards_symbol,
          "|rank=", genes$genecards_rank,
          "|score=", .format_score(genes$genecards_relevance_score)
        )
      } else character()
    ),
    max_relevance_score_among_overlapping_genes = maximum_score,
    highest_relevance_gene_symbols = .collapse(
      if (contains) genes$genecards_symbol[maximum_rows] else character()
    ),
    highest_relevance_gene_ranks = .collapse(
      if (contains) genes$genecards_rank[maximum_rows] else character()
    ),
    max_relevance_score_among_fully_contained_genes =
      fully_contained_maximum_score,
    highest_relevance_fully_contained_gene_symbols = .collapse(
      genes$genecards_symbol[fully_contained_maximum_rows]
    ),
    coordinate_system = "GRCh37/hg19",
    overlap_rule = "closed_gene_body_overlap",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = NULL
  )
}

region_summary <- do.call(rbind, summary_parts)
rownames(region_summary) <- NULL
score_levels <- sort(unique(
  region_summary$max_relevance_score_among_overlapping_genes[
    is.finite(region_summary$max_relevance_score_among_overlapping_genes)
  ]
), decreasing = TRUE)
ranked_regions <- region_summary
ranked_regions$relevance_rank <- match(
  ranked_regions$max_relevance_score_among_overlapping_genes,
  score_levels
)
ranked_regions$ranking_status <- ifelse(
  ranked_regions$contains_genecards_top1000_gene,
  "ranked_by_maximum_overlapping_gene_score",
  "not_ranked_no_top1000_gene"
)
ranked_regions <- ranked_regions[order(
  is.na(ranked_regions$relevance_rank),
  ranked_regions$relevance_rank,
  match(ranked_regions$all_10_class, class_levels),
  ranked_regions$chromosome,
  ranked_regions$start
), ]
ranked_regions <- ranked_regions[, c(
  "relevance_rank", "ranking_status",
  setdiff(names(ranked_regions), c("relevance_rank", "ranking_status"))
)]
rownames(ranked_regions) <- NULL
long <- if (length(long_parts)) {
  do.call(rbind, long_parts)
} else {
  stop("No GeneCards top-1000 overlaps were found.", call. = FALSE)
}
rownames(long) <- NULL

pair_key <- paste(long$region_id, long$genecards_rank, sep = "|")
if (anyDuplicated(pair_key)) {
  stop(
    "A GeneCards input maps to multiple overlapping loci in one region; ",
    "the region-gene definition must be revisited.",
    call. = FALSE
  )
}

top500_long <- read.delim(
  top500_long_file, check.names = FALSE, stringsAsFactors = FALSE
)
regression_columns <- setdiff(names(long), "fraction_gene_span_overlapped")
top500_long_regression <-
  identical(names(top500_long), names(long)) &&
  identical(
    .normalize_table_for_regression(
      long[long$genecards_rank <= 500L, regression_columns, drop = FALSE]
    ),
    .normalize_table_for_regression(
      top500_long[, regression_columns, drop = FALSE]
    )
  )

.summarize_class <- function(class_name) {
  if (identical(class_name, "All_classes")) {
    region_subset <- region_summary
    hit_subset <- long
    class_label <- "All classes"
  } else {
    region_subset <- region_summary[
      region_summary$all_10_class == class_name, , drop = FALSE
    ]
    hit_subset <- long[long$all_10_class == class_name, , drop = FALSE]
    class_label <- unname(display_class[class_name])
  }

  unique_hit_rows <- hit_subset[!duplicated(hit_subset$genecards_rank), ]
  unique_hit_rows <- unique_hit_rows[order(
    -unique_hit_rows$relevance_score,
    unique_hit_rows$genecards_rank,
    unique_hit_rows$genecards_symbol
  ), ]
  class_max <- if (nrow(hit_subset)) max(hit_subset$relevance_score) else NA_real_
  overlap_region_scores <-
    region_subset$max_relevance_score_among_overlapping_genes
  overlap_region_scores <- overlap_region_scores[is.finite(overlap_region_scores)]
  strict_region_scores <-
    region_subset$max_relevance_score_among_fully_contained_genes
  strict_region_scores <- strict_region_scores[is.finite(strict_region_scores)]

  data.frame(
    all_10_class = class_name,
    display_class = class_label,
    total_regions = nrow(region_subset),
    regions_containing_top1000_gene = sum(
      region_subset$contains_genecards_top1000_gene
    ),
    regions_not_containing_top1000_gene = sum(
      !region_subset$contains_genecards_top1000_gene
    ),
    fraction_regions_containing_top1000_gene = mean(
      region_subset$contains_genecards_top1000_gene
    ),
    region_gene_overlap_pairs = nrow(hit_subset),
    unique_genecards_top1000_genes = nrow(unique_hit_rows),
    unique_top1000_gene_symbols_by_relevance = .collapse(
      unique_hit_rows$genecards_symbol
    ),
    mean_region_max_score_among_overlapping_genes =
      mean(overlap_region_scores),
    median_region_max_score_among_overlapping_genes =
      median(overlap_region_scores),
    regions_with_fully_contained_top1000_gene =
      length(strict_region_scores),
    mean_region_max_score_among_fully_contained_genes =
      mean(strict_region_scores),
    median_region_max_score_among_fully_contained_genes =
      median(strict_region_scores),
    class_max_relevance_score_among_overlapping_genes = class_max,
    class_highest_relevance_gene_symbols = .collapse(
      if (nrow(hit_subset)) {
        unique(hit_subset$genecards_symbol[hit_subset$relevance_score == class_max])
      } else character()
    ),
    stringsAsFactors = FALSE
  )
}

class_summary <- do.call(
  rbind, lapply(c(class_levels, "All_classes"), .summarize_class)
)
rownames(class_summary) <- NULL

summary_hit_counts <- setNames(
  class_summary$regions_containing_top1000_gene[match(
    class_levels, class_summary$all_10_class
  )],
  class_levels
)

output_checks <- c(
  one_summary_row_per_region = nrow(region_summary) == 27L &&
    !anyDuplicated(region_summary$region_id),
  every_source_region_retained = setequal(
    region_summary$region_id, regions$region_id
  ),
  long_pairs_unique = !anyDuplicated(pair_key),
  top500_long_regression = top500_long_regression,
  ranking_has_every_region_once = nrow(ranked_regions) == 27L &&
    !anyDuplicated(ranked_regions$region_id),
  ranking_scores_descend = all(diff(
    ranked_regions$max_relevance_score_among_overlapping_genes[
      is.finite(ranked_regions$max_relevance_score_among_overlapping_genes)
    ]
  ) <= 0),
  no_hit_regions_are_unranked = all(is.na(
    ranked_regions$relevance_rank[
      !ranked_regions$contains_genecards_top1000_gene
    ]
  )),
  long_intervals_really_overlap = all(
    long$gene_start_1based <= long$end &
      long$gene_end_1based >= long$start &
      long$overlap_bp >= 1L
  ),
  long_contains_no_functional_elements = !any(
    long$genecards_type == "Functional Element"
  ),
  region_counts_match_long = all(
    region_summary$n_genecards_top1000_genes ==
      as.integer(table(factor(
        long$region_id,
        levels = region_summary$region_id
      )))
  ),
  fully_contained_counts_match_long = all(
    region_summary$n_fully_contained_top1000_genes ==
      as.integer(table(factor(
        long$region_id[long$fully_contained_in_region],
        levels = region_summary$region_id
      )))
  ),
  no_hit_regions_have_no_maximum = all(
    is.na(region_summary$max_relevance_score_among_overlapping_genes[
      !region_summary$contains_genecards_top1000_gene
    ])
  ),
  hit_regions_have_correct_maximum = all(vapply(
    region_summary$region_id[region_summary$contains_genecards_top1000_gene],
    function(region_id) {
      observed <- region_summary$max_relevance_score_among_overlapping_genes[
        region_summary$region_id == region_id
      ]
      expected <- max(long$relevance_score[long$region_id == region_id])
      identical(as.numeric(observed), as.numeric(expected))
    },
    logical(1L)
  )),
  no_fully_contained_gene_means_no_strict_maximum = all(
    is.na(region_summary$max_relevance_score_among_fully_contained_genes[
      region_summary$n_fully_contained_top1000_genes == 0L
    ])
  ),
  fully_contained_maximum_is_correct = all(vapply(
    region_summary$region_id[
      region_summary$n_fully_contained_top1000_genes > 0L
    ],
    function(region_id) {
      observed <-
        region_summary$max_relevance_score_among_fully_contained_genes[
          region_summary$region_id == region_id
        ]
      expected <- max(long$relevance_score[
        long$region_id == region_id & long$fully_contained_in_region
      ])
      identical(as.numeric(observed), as.numeric(expected))
    },
    logical(1L)
  )),
  class_region_totals = identical(
    class_summary$total_regions[match(class_levels, class_summary$all_10_class)],
    unname(as.integer(expected_class_counts))
  ),
  all_class_total_27 = class_summary$total_regions[
    class_summary$all_10_class == "All_classes"
  ] == 27L,
  class_hit_totals_cover_all_hits = sum(summary_hit_counts) ==
    sum(region_summary$contains_genecards_top1000_gene)
)
if (!all(output_checks)) {
  stop(
    "Output validation failed: ",
    paste(names(output_checks)[!output_checks], collapse = ", "),
    call. = FALSE
  )
}

.write_tsv(region_summary, summary_file)
.write_tsv(long, long_file)
.write_tsv(class_summary, class_summary_file)
.write_tsv(ranked_regions, ranking_file)

all_checks <- c(region_checks, gene_checks, output_checks)
metadata <- c(
  "analysis=AFR all-10-seed stable Height regions versus GeneCards Height top-1000 genes",
  "mapping_snapshot_date=2026-08-30",
  "target_ancestry=AFR",
  "knockoff_copies_M=1",
  "seeds=100,200,300,400,500,600,700,800,900,1000",
  "q_threshold=0.1",
  "region_set=union of GhostKnockoff and AnnoGK regions detected in all 10 seeds",
  paste0("implementation=", script_file),
  paste0("implementation_md5=", .md5(script_file)),
  paste0("region_file=", region_file),
  paste0("region_file_md5=", .md5(region_file)),
  paste0("gene_file=", gene_file),
  paste0("gene_file_md5=", .md5(gene_file)),
  paste0("gene_mapping_marker=", gene_mapping_marker),
  paste0("gene_mapping_marker_md5=", .md5(gene_mapping_marker)),
  paste0("top500_overlap_long_file=", top500_long_file),
  paste0("top500_overlap_long_file_md5=", .md5(top500_long_file)),
  paste0("top500_overlap_marker=", top500_overlap_marker),
  paste0("top500_overlap_marker_md5=", .md5(top500_overlap_marker)),
  "coordinate_system=GRCh37/hg19",
  "interval_convention=1-based closed",
  paste(
    "overlap_rule=a GeneCards gene matches a region when gene_start <=",
    "region_end and gene_end >= region_start"
  ),
  paste(
    "containment_note=fully_contained_in_region is also reported, but region",
    "membership is defined by any closed gene-body overlap"
  ),
  paste(
    "score_rule=per-region maximum is the maximum GeneCards Relevance Score",
    "among overlapping mapped gene/RNA inputs; it is missing for regions",
    "with no top-1000 overlap"
  ),
  paste(
    "ranking_rule=dense rank of per-region maximum overlapping-gene score;",
    "ties share a rank and no-hit regions are not assigned zero or a rank"
  ),
  "functional_element_rule=the ten GeneCards functional elements are excluded from gene overlap",
  "long_table_rule=one row per unique region and GeneCards ranked gene input",
  paste0(
    "validation=",
    paste(names(all_checks), all_checks, sep = ":", collapse = ";")
  )
)
.write_lines(metadata, metadata_file)

output_files <- c(
  summary_file, long_file, class_summary_file, ranking_file, metadata_file
)
completion <- c(
  "status=complete",
  "scope=AFR_all10_seed_regions",
  "coordinate_system=GRCh37/hg19",
  paste0("regions=", nrow(region_summary)),
  paste0("regions_with_top1000_gene=", sum(region_summary$contains_genecards_top1000_gene)),
  paste0("regions_without_top1000_gene=", sum(!region_summary$contains_genecards_top1000_gene)),
  paste0("region_gene_pairs=", nrow(long)),
  paste0("unique_top1000_genes=", length(unique(long$genecards_rank))),
  paste0("validation_checks=", length(all_checks)),
  "all_validations_passed=true",
  paste0("implementation_md5=", .md5(script_file)),
  paste0(
    "output_md5=",
    paste(
      basename(output_files),
      vapply(output_files, .md5, character(1L)),
      sep = ":", collapse = ";"
    )
  )
)
.write_lines(completion, completion_file)

cat("Completed AFR all-10 region overlap with GeneCards Height top-1000 genes.\n")
print(class_summary[, c(
  "display_class", "total_regions", "regions_containing_top1000_gene",
  "regions_not_containing_top1000_gene", "unique_genecards_top1000_genes",
  "mean_region_max_score_among_overlapping_genes",
  "median_region_max_score_among_overlapping_genes",
  "class_max_relevance_score_among_overlapping_genes",
  "class_highest_relevance_gene_symbols"
)], row.names = FALSE)
