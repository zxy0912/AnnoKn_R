if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
setwd(repro_path("gtex_validation_workspace"))
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1L]] else {
  "gtex_method_specific_eGene_gtex_v8_comparison_q01_old_seed_PET"
}

overall_path <- file.path(input_dir, "gtex_v8_overlap_by_method.tsv")
tissue_path <- file.path(
  input_dir,
  "gtex_v8_overlap_by_method_and_tissue.tsv"
)

overall <- read.delim(overall_path, check.names = FALSE)
by_tissue <- read.delim(tissue_path, check.names = FALSE)

required_groups <- c("AnnoKn-only", "Vanilla-Knockoff-only")
if (!setequal(overall$method_specific_group, required_groups)) {
  stop("Expected exactly the AnnoKn-only and Vanilla-Knockoff-only groups")
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator
  c(lower = center - half_width, upper = center + half_width)
}

intervals <- t(mapply(
  wilson_interval,
  overall$gtex_v8_standard_egene,
  overall$n_pairs
))
overall$support_percent <- 100 * overall$support_fraction
overall$wilson_95_lower <- intervals[, "lower"]
overall$wilson_95_upper <- intervals[, "upper"]

ann <- overall[overall$method_specific_group == "AnnoKn-only", ]
van <- overall[overall$method_specific_group == "Vanilla-Knockoff-only", ]

contingency <- matrix(
  c(
    ann$gtex_v8_standard_egene,
    ann$not_gtex_v8_standard_egene,
    van$gtex_v8_standard_egene,
    van$not_gtex_v8_standard_egene
  ),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(
    group = c("AnnoKn-only", "Vanilla-Knockoff-only"),
    gtex_v8 = c("eGene", "not_eGene")
  )
)

prop_result <- prop.test(
  x = c(ann$gtex_v8_standard_egene, van$gtex_v8_standard_egene),
  n = c(ann$n_pairs, van$n_pairs),
  correct = FALSE
)
fisher_result <- fisher.test(contingency)

rr <- ann$support_fraction / van$support_fraction
rr_se <- sqrt(
  1 / ann$gtex_v8_standard_egene - 1 / ann$n_pairs +
    1 / van$gtex_v8_standard_egene - 1 / van$n_pairs
)
rr_ci <- exp(log(rr) + c(-1, 1) * qnorm(0.975) * rr_se)

comparison <- data.frame(
  contrast = "AnnoKn-only minus Vanilla-Knockoff-only",
  annokn_only_supported = ann$gtex_v8_standard_egene,
  annokn_only_total = ann$n_pairs,
  annokn_only_fraction = ann$support_fraction,
  vanilla_only_supported = van$gtex_v8_standard_egene,
  vanilla_only_total = van$n_pairs,
  vanilla_only_fraction = van$support_fraction,
  fraction_difference = ann$support_fraction - van$support_fraction,
  difference_95_lower = prop_result$conf.int[[1L]],
  difference_95_upper = prop_result$conf.int[[2L]],
  two_sample_proportion_p = prop_result$p.value,
  risk_ratio = rr,
  risk_ratio_95_lower = rr_ci[[1L]],
  risk_ratio_95_upper = rr_ci[[2L]],
  odds_ratio = unname(fisher_result$estimate),
  odds_ratio_95_lower = fisher_result$conf.int[[1L]],
  odds_ratio_95_upper = fisher_result$conf.int[[2L]],
  fisher_exact_p = fisher_result$p.value,
  stringsAsFactors = FALSE
)

tissue_wide <- reshape(
  by_tissue[, c(
    "tissue",
    "method_specific_group",
    "n_pairs",
    "gtex_v8_standard_egene",
    "support_fraction"
  )],
  idvar = "tissue",
  timevar = "method_specific_group",
  direction = "wide"
)
names(tissue_wide) <- gsub("[.]AnnoKn-only$", "_annokn_only", names(tissue_wide))
names(tissue_wide) <- gsub(
  "[.]Vanilla-Knockoff-only$",
  "_vanilla_knockoff_only",
  names(tissue_wide)
)
tissue_wide$fraction_difference <- with(
  tissue_wide,
  support_fraction_annokn_only - support_fraction_vanilla_knockoff_only
)
tissue_wide <- tissue_wide[order(tissue_wide$tissue), ]

stratified_array <- array(
  0,
  dim = c(2L, 2L, nrow(tissue_wide)),
  dimnames = list(
    group = c("AnnoKn-only", "Vanilla-Knockoff-only"),
    gtex_v8 = c("eGene", "not_eGene"),
    tissue = tissue_wide$tissue
  )
)
stratified_array["AnnoKn-only", "eGene", ] <-
  tissue_wide$gtex_v8_standard_egene_annokn_only
stratified_array["AnnoKn-only", "not_eGene", ] <-
  tissue_wide$n_pairs_annokn_only -
    tissue_wide$gtex_v8_standard_egene_annokn_only
stratified_array["Vanilla-Knockoff-only", "eGene", ] <-
  tissue_wide$gtex_v8_standard_egene_vanilla_knockoff_only
stratified_array["Vanilla-Knockoff-only", "not_eGene", ] <-
  tissue_wide$n_pairs_vanilla_knockoff_only -
    tissue_wide$gtex_v8_standard_egene_vanilla_knockoff_only

mh_result <- mantelhaen.test(stratified_array, correct = FALSE)
combined_tissue_weights <- with(
  tissue_wide,
  (n_pairs_annokn_only + n_pairs_vanilla_knockoff_only) /
    sum(n_pairs_annokn_only + n_pairs_vanilla_knockoff_only)
)
standardized_annokn <- sum(
  combined_tissue_weights * tissue_wide$support_fraction_annokn_only
)
standardized_vanilla <- sum(
  combined_tissue_weights *
    tissue_wide$support_fraction_vanilla_knockoff_only
)
stratified_comparison <- data.frame(
  strata = "11 GTEx tissues",
  common_odds_ratio_annokn_vs_vanilla = unname(mh_result$estimate),
  common_odds_ratio_95_lower = mh_result$conf.int[[1L]],
  common_odds_ratio_95_upper = mh_result$conf.int[[2L]],
  mantel_haenszel_p = mh_result$p.value,
  combined_weight_standardized_annokn_fraction = standardized_annokn,
  combined_weight_standardized_vanilla_fraction = standardized_vanilla,
  standardized_fraction_difference_annokn_minus_vanilla =
    standardized_annokn - standardized_vanilla,
  equal_tissue_mean_annokn_fraction = mean(
    tissue_wide$support_fraction_annokn_only
  ),
  equal_tissue_mean_vanilla_fraction = mean(
    tissue_wide$support_fraction_vanilla_knockoff_only
  ),
  equal_tissue_mean_difference_annokn_minus_vanilla = mean(
    tissue_wide$fraction_difference
  ),
  tissues_annokn_higher = sum(tissue_wide$fraction_difference > 0),
  tissues_tied = sum(tissue_wide$fraction_difference == 0),
  tissues_vanilla_higher = sum(tissue_wide$fraction_difference < 0),
  stringsAsFactors = FALSE
)

write.table(
  overall,
  file.path(input_dir, "gtex_v8_overlap_by_method_with_ci.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  comparison,
  file.path(input_dir, "gtex_v8_overlap_global_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  tissue_wide,
  file.path(input_dir, "gtex_v8_overlap_tissue_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  stratified_comparison,
  file.path(input_dir, "gtex_v8_overlap_tissue_stratified_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

print(overall, row.names = FALSE)
print(comparison, row.names = FALSE)
print(stratified_comparison, row.names = FALSE)
