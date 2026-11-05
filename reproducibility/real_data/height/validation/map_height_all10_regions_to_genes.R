if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
#!/usr/bin/env Rscript

# Map the union of AFR regions detected in all ten M=1 seeds by either
# GhostKnockoff or AnnoGK to GENCODE v19 genes on GRCh37.
#
# Mapping is intentionally descriptive:
#   1. a gene-body overlap means [gene_start, gene_end] intersects the region;
#   2. protein-coding genes are summarized separately from all GENCODE genes;
#   3. genes with a TSS inside the region are also reported;
#   4. no mapped gene is labelled causal.

options(stringsAsFactors = FALSE)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- repro_path("height_validation_workspace")

input_file <- file.path(
  script_dir,
  "height_ten_seed_region_stability_q0p1",
  "AFR_all10_union_regions.tsv"
)
gtf_file <- file.path(script_dir, "reference", "gencode.v19.annotation.gtf.gz")
output_dir <- file.path(script_dir, "height_ten_seed_region_stability_q0p1")

summary_file <- file.path(
  output_dir,
  "AFR_all10_union_regions_gene_mapping_gencode_v19.tsv"
)
all_long_file <- file.path(
  output_dir,
  "AFR_all10_union_regions_gene_mapping_gencode_v19_all_genes_long.tsv"
)
pc_long_file <- file.path(
  output_dir,
  "AFR_all10_union_regions_gene_mapping_gencode_v19_protein_coding_long.tsv"
)
metadata_file <- file.path(
  output_dir,
  "AFR_all10_union_regions_gene_mapping_gencode_v19_metadata.txt"
)

stopifnot(file.exists(input_file), file.exists(gtf_file))

regions <- read.delim(input_file, check.names = FALSE)
required_region_columns <- c(
  "region_id", "ancestry", "chromosome", "start", "end", "all_10_class"
)
stopifnot(all(required_region_columns %in% names(regions)))
stopifnot(nrow(regions) == 27L)
stopifnot(all(regions$start <= regions$end))

regions$chromosome <- as.character(regions$chromosome)
regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)
regions$chromosome_gtf <- paste0("chr", regions$chromosome)

class_levels <- c("Both", "AnnoGK_only", "GhostKnockoff_only")
stopifnot(all(regions$all_10_class %in% class_levels))
regions$class_order <- match(regions$all_10_class, class_levels)
regions <- regions[order(regions$class_order, as.numeric(regions$chromosome), regions$start), ]
row.names(regions) <- NULL

# Read only GTF gene records in bounded chunks so the full transcript-level GTF
# never has to be held in memory.
con <- gzfile(gtf_file, open = "rt")
on.exit(close(con), add = TRUE)
gene_lines <- list()
chunk_index <- 0L
repeat {
  chunk <- readLines(con, n = 100000L, warn = FALSE)
  if (!length(chunk)) break
  keep <- grepl("\tgene\t", chunk, fixed = TRUE)
  if (any(keep)) {
    chunk_index <- chunk_index + 1L
    gene_lines[[chunk_index]] <- chunk[keep]
  }
}
close(con)
on.exit(NULL, add = FALSE)
gene_lines <- unlist(gene_lines, use.names = FALSE)
stopifnot(length(gene_lines) > 50000L)

fields <- strsplit(gene_lines, "\t", fixed = TRUE)
stopifnot(all(lengths(fields) == 9L))
gtf <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
names(gtf) <- c(
  "chromosome", "source", "feature", "gene_start", "gene_end",
  "score", "strand", "frame", "attributes"
)
gtf$gene_start <- as.numeric(gtf$gene_start)
gtf$gene_end <- as.numeric(gtf$gene_end)

extract_gtf_attribute <- function(x, key) {
  pattern <- paste0('.*(?:^|; )', key, ' "([^"]+)".*')
  found <- grepl(pattern, x, perl = TRUE)
  answer <- rep(NA_character_, length(x))
  answer[found] <- sub(pattern, "\\1", x[found], perl = TRUE)
  answer
}

gtf$gene_id <- extract_gtf_attribute(gtf$attributes, "gene_id")
gtf$gene_id_base <- sub("\\..*$", "", gtf$gene_id)
gtf$gene_name <- extract_gtf_attribute(gtf$attributes, "gene_name")
gtf$gene_type <- extract_gtf_attribute(gtf$attributes, "gene_type")
gtf$tss <- ifelse(gtf$strand == "+", gtf$gene_start, gtf$gene_end)

stopifnot(!anyNA(gtf$gene_id), !anyNA(gtf$gene_name), !anyNA(gtf$gene_type))
stopifnot(!anyDuplicated(gtf$gene_id))

region_columns <- c(
  "region_id", "ancestry", "chromosome", "start", "end", "all_10_class",
  "ghost_n_seeds", "annogk_n_seeds"
)
gene_columns <- c(
  "gene_id", "gene_id_base", "gene_name", "gene_type", "gene_start",
  "gene_end", "strand", "tss"
)

all_long <- vector("list", nrow(regions))
pc_long <- vector("list", nrow(regions))
summary_rows <- vector("list", nrow(regions))

collapse_sorted <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x)) paste(x, collapse = ";") else ""
}

for (i in seq_len(nrow(regions))) {
  region <- regions[i, ]
  same_chr <- gtf$chromosome == region$chromosome_gtf
  overlaps <- same_chr & gtf$gene_start <= region$end & gtf$gene_end >= region$start
  mapped <- gtf[overlaps, gene_columns, drop = FALSE]
  mapped <- mapped[order(mapped$gene_start, mapped$gene_end, mapped$gene_name), ]

  all_long[[i]] <- cbind(
    region[rep(1L, nrow(mapped)), region_columns, drop = FALSE],
    overlap_type = rep("gene_body_overlap", nrow(mapped)),
    mapped,
    row.names = NULL
  )

  mapped_pc <- mapped[mapped$gene_type == "protein_coding", , drop = FALSE]
  pc_long[[i]] <- cbind(
    region[rep(1L, nrow(mapped_pc)), region_columns, drop = FALSE],
    overlap_type = rep("gene_body_overlap", nrow(mapped_pc)),
    mapped_pc,
    row.names = NULL
  )

  tss_inside_all <- mapped$tss >= region$start & mapped$tss <= region$end
  tss_inside_pc <- mapped_pc$tss >= region$start & mapped_pc$tss <= region$end

  # This fallback is useful only if a region has no protein-coding TSS inside it.
  chr_pc <- gtf[same_chr & gtf$gene_type == "protein_coding", , drop = FALSE]
  tss_distance <- ifelse(
    chr_pc$tss < region$start,
    region$start - chr_pc$tss,
    ifelse(chr_pc$tss > region$end, chr_pc$tss - region$end, 0)
  )
  min_tss_distance <- min(tss_distance)
  nearest_pc <- chr_pc$gene_name[tss_distance == min_tss_distance]

  summary_rows[[i]] <- data.frame(
    region_id = region$region_id,
    ancestry = region$ancestry,
    chromosome = region$chromosome,
    start = region$start,
    end = region$end,
    all_10_class = region$all_10_class,
    ghost_n_seeds = region$ghost_n_seeds,
    annogk_n_seeds = region$annogk_n_seeds,
    n_all_gene_body_overlaps = nrow(mapped),
    n_protein_coding_gene_body_overlaps = nrow(mapped_pc),
    protein_coding_gene_body_overlaps = collapse_sorted(mapped_pc$gene_name),
    n_all_gene_tss_inside = sum(tss_inside_all),
    n_protein_coding_gene_tss_inside = sum(tss_inside_pc),
    protein_coding_gene_tss_inside = collapse_sorted(mapped_pc$gene_name[tss_inside_pc]),
    nearest_protein_coding_tss_gene = collapse_sorted(nearest_pc),
    nearest_protein_coding_tss_distance_bp = min_tss_distance,
    stringsAsFactors = FALSE
  )
}

all_long <- do.call(rbind, all_long)
pc_long <- do.call(rbind, pc_long)
summary_table <- do.call(rbind, summary_rows)

stopifnot(nrow(summary_table) == 27L)
stopifnot(all(summary_table$n_all_gene_body_overlaps > 0L))
stopifnot(all(summary_table$n_protein_coding_gene_body_overlaps > 0L))
stopifnot(all(summary_table$n_protein_coding_gene_tss_inside > 0L))
stopifnot(sum(summary_table$all_10_class == "Both") == 9L)
stopifnot(sum(summary_table$all_10_class == "AnnoGK_only") == 12L)
stopifnot(sum(summary_table$all_10_class == "GhostKnockoff_only") == 6L)

write.table(summary_table, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_long, all_long_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pc_long, pc_long_file, sep = "\t", quote = FALSE, row.names = FALSE)

metadata <- c(
  "Analysis: AFR all-10-seed Height-region to gene mapping",
  "Genome build: GRCh37 / hg19",
  "Gene annotation: GENCODE release 19 comprehensive gene annotation, reference chromosomes",
  "GENCODE assembly: GRCh37.p13",
  "GENCODE release date: 2013-12-05",
  "GENCODE URL: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
  "GENCODE SHA256: 9909e8123a40eca090c3f8e501633f7cc65cd6c4adb83b95009745eca3564d70",
  paste0("Input regions: ", normalizePath(input_file)),
  paste0("Input annotation: ", normalizePath(gtf_file)),
  "Region coordinates and GTF gene coordinates are treated as 1-based closed intervals.",
  "Gene-body overlap rule: gene_start <= region_end AND gene_end >= region_start.",
  "TSS rule: gene_start on the + strand and gene_end on the - strand.",
  "The protein-coding summary uses GENCODE gene_type == protein_coding.",
  "A mapped/nearest gene is a positional annotation, not a causal-gene assignment.",
  paste0("Number of union regions: ", nrow(summary_table)),
  paste0("All GENCODE gene-body overlap rows: ", nrow(all_long)),
  paste0("Protein-coding gene-body overlap rows: ", nrow(pc_long))
)
writeLines(metadata, metadata_file)

message("Wrote: ", summary_file)
message("Wrote: ", all_long_file)
message("Wrote: ", pc_long_file)
message("Wrote: ", metadata_file)
