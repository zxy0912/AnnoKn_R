if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
# Reproduce the supplied four-method SCZ region analysis with coordinate IDs.
# Usage: Rscript SCZ/summarize_scz_seed_regions.R [result_dir] [output_dir]
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[1L] else
  repro_path("scz_validation_workspace", "result")
output_dir <- if (length(args) >= 2L) args[2L] else
  file.path(dirname(input_dir), "hundred_seed_region_stability_q0p1")
seeds <- 100L * seq_len(100L)
chromosomes <- setdiff(seq_len(22L), 21L)
ancestries <- c("EUR", "EAS")
methods <- c("GhostKnockoff", "AnnoGK", "Transfer", "AdaKn")
threshold <- 0.1
manifest <- expand.grid(seed = seeds, chr = chromosomes, ancestry = ancestries,
                        stringsAsFactors = FALSE)
other <- ifelse(manifest$ancestry == "EUR", "EAS_AFR_LAT", "EUR_AFR_LAT")
manifest$file <- file.path(input_dir, sprintf(
  "%s_result_adakn_10_median_r2_%s_chr_%d_M_1_%d.RData",
  manifest$ancestry, other, manifest$chr, manifest$seed))
stopifnot(nrow(manifest) == 4200L)
missing <- manifest$file[!file.exists(manifest$file)]
if (length(missing)) stop("Missing files: ", paste(missing, collapse = "\n"))
records <- vector("list", nrow(manifest))
reference_regions <- list()
for (j in seq_len(nrow(manifest))) {
  row <- manifest[j, ]
  env <- new.env(parent = baseenv())
  load(row$file, envir = env)
  stopifnot(exists("result", envir = env, inherits = FALSE))
  result <- env$result
  fields <- c("risk_region", "snp", "q_gk", "q_gkanno", "q_transfer", "rej_adakn")
  stopifnot(all(fields %in% names(result)), result$chr == row$chr)
  n <- length(result$risk_region)
  stopifnot(n > 0L, all(lengths(result[fields]) == n))
  coords <- t(vapply(result$risk_region, function(x) {
    z <- as.numeric(unlist(x, use.names = FALSE))
    stopifnot(length(z) == 2L, all(is.finite(z)), z[1] <= z[2])
    z
  }, numeric(2)))
  ids <- paste(row$chr, coords[, 1], coords[, 2], sep = "_")
  stopifnot(!anyDuplicated(ids))
  key <- paste(row$ancestry, row$chr, sep = "_")
  if (is.null(reference_regions[[key]])) reference_regions[[key]] <- ids
  stopifnot(identical(reference_regions[[key]], ids))
  hits <- matrix(FALSE, nrow = n, ncol = 4L, dimnames = list(NULL, methods))
  for (i in seq_len(n)) {
    nsnp <- nrow(result$snp[[i]])
    for (k in seq_len(3L)) {
      q <- result[[c("q_gk", "q_gkanno", "q_transfer")[k]]][[i]]
      stopifnot(is.numeric(q), length(q) == nsnp, !anyNA(q))
      hits[i, k] <- any(q <= threshold)
    }
    d <- result$rej_adakn[[i]]
    # Preserve the user's stored AdaKn rejection-index rule; do not refit it.
    if (length(d)) stopifnot(is.numeric(d), !anyNA(d),
                            all(d == as.integer(d)), all(d >= 1L & d <= nsnp))
    hits[i, 4L] <- length(d) > 0L
  }
  records[[j]] <- data.frame(ancestry = row$ancestry, seed = row$seed,
    chr = row$chr, region_index = seq_len(n), start = coords[, 1],
    end = coords[, 2], region_id = ids, hits, check.names = FALSE)
  if (j %% 500L == 0L) message("Read ", j, "/", nrow(manifest), " files")
}
detections <- do.call(rbind, records)
stopifnot(!anyDuplicated(detections[c("ancestry", "seed", "region_id")]))
group_columns <- c("ancestry", "chr", "region_index", "start", "end", "region_id")
recurrence <- aggregate(detections[methods], detections[group_columns], sum)
observed <- aggregate(list(n_seeds_observed = detections$seed),
                      detections[group_columns], length)
stopifnot(all(observed$n_seeds_observed == length(seeds)))
recurrence <- merge(recurrence, observed, by = group_columns, sort = FALSE)
recurrence <- recurrence[order(recurrence$ancestry, recurrence$chr, recurrence$start), ]
summary <- do.call(rbind, lapply(methods, function(method) {
  sets <- lapply(ancestries, function(a) {
    x <- recurrence[recurrence$ancestry == a, ]
    list(union = x$region_id[x[[method]] >= 1L],
         intersection = x$region_id[x[[method]] == length(seeds)])
  })
  data.frame(method = method,
    EUR_union = length(sets[[1]]$union), EAS_union = length(sets[[2]]$union),
    common_union = length(intersect(sets[[1]]$union, sets[[2]]$union)),
    EUR_all100 = length(sets[[1]]$intersection),
    EAS_all100 = length(sets[[2]]$intersection),
    common_all100 = length(intersect(sets[[1]]$intersection, sets[[2]]$intersection)))
}))
# Check the original chr_i/table implementation independently against coordinate counts.
for (a in ancestries) for (method in methods) {
  x <- detections[detections$ancestry == a & detections[[method]], ]
  original_counts <- table(paste(x$chr, x$region_index, sep = "_"))
  r <- recurrence[recurrence$ancestry == a, ]
  stopifnot(length(original_counts) == sum(r[[method]] >= 1L),
            sum(original_counts == 100L) == sum(r[[method]] == 100L))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(output_dir, name),
  sep = "\t", quote = FALSE, row.names = FALSE)
write_tsv(summary, "method_region_summary.tsv")
write_tsv(recurrence, "region_seed_recurrence.tsv")
write_tsv(detections, "region_by_seed_detections.tsv")
stable <- do.call(rbind, lapply(methods, function(method) {
  x <- recurrence[recurrence[[method]] == length(seeds), ]
  x$method <- rep(method, nrow(x))
  x
}))
write_tsv(stable, "all100_regions.tsv")
manifest$md5 <- unname(tools::md5sum(manifest$file))
write_tsv(manifest, "input_manifest.tsv")
capture.output(list(input_dir = normalizePath(input_dir), seeds = seeds,
  chromosomes = chromosomes, M = 1L, q_threshold = threshold,
  AdaKn_rule = "length(result$rej_adakn[[i]]) > 0; stored rejection indices",
  region_definition = "Exact chromosome/start/end; any selected SNP per region per seed",
  coordinate_order_consistent_across_seeds = TRUE,
  all_regions_observed_in_100_seeds = TRUE,
  original_chr_index_counts_match = TRUE,
  session = sessionInfo()), file = file.path(output_dir, "metadata.txt"))
print(summary, row.names = FALSE)
message("Outputs: ", normalizePath(output_dir))
