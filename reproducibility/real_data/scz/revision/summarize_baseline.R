if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)


############## check if multiple knockoff copies


find1 <- c()
find1_anno <- c()
find2 <- c()
find2_anno <- c()

number_gk_eur <- c()
number_annogk_eur <- c()
number_gk_eas <- c()
number_annogk_eas <- c()

#

all_seed = 100*(1:100)
# c(1,12,10,10000,123456, 100, 1000, 123, 1234, 12345)
for(seed in all_seed){
  region_1 <- numeric()
  region_anno_1 <- numeric()
  region_2 <- numeric()
  region_anno_2 <- numeric()
  lambdas1 <- list()
  lambdas2 <- list()

  threshold = 0.1
  sum1 = 0
  sum2 = 0
  chrid = '1'
  M = 1
  # seed = '12345' for dss


  chrlist = (1:22)[c(-21)]
  for(chrid in chrlist){
    print(paste0("chrid", chrid))
    ancestry0 = 'EUR'
    ancestry1 = c('EAS',"AFR","LAT")
    # ancestry1 = c('EAS')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("scz_baseline_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)

    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)

    for(i in 1:n_region){
      # print(i)
      snp = result$snp[[i]]
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
        find1 <- append(find1, paste(chrid, i, sep = '_'))
        number_gk_eur <- append(number_gk_eur, snp$id[which(a <= threshold)])
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
        find1_anno <- append(find1_anno, paste(chrid, i, sep = '_'))
        number_annogk_eur <- append(number_annogk_eur, snp$id[which(b <= threshold)])
      }
    }


    ancestry0 = 'EAS'
    ancestry1 = c('EUR',"AFR","LAT")
    # ancestry1 = c('EUR')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("scz_baseline_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)
    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)


    for(i in 1:n_region){
      # print(i)
      snp = result$snp[[i]]
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
        find2 <- append(find2, paste(chrid, i, sep = '_'))
        number_gk_eas <- append(number_gk_eas, snp$id[which(a <= threshold)])
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
        find2_anno <- append(find2_anno, paste(chrid, i, sep = '_'))
        number_annogk_eas <- append(number_annogk_eas, snp$id[which(b <= threshold)])
      }
    }

  }

  nrow(region_1)
  nrow(region_2)

  nrow(region_anno_1)
  nrow(region_anno_2)
}


length(unique(number_gk_eur))
length(unique(number_annogk_eur))
length(unique(number_gk_eas))
length(unique(number_annogk_eas))

table(table(number_gk_eur))
table(table(number_annogk_eur))
table(table(number_gk_eas))
table(table(number_annogk_eas))



repli = 4

sum(table(find1) > repli)
sum(table(find1_anno) > repli)

sum(table(find2) > repli)
sum(table(find2_anno) > repli)


length(table(find1))
length(table(find1_anno))
length(table(find2))
length(table(find2_anno))



length(intersect(unique(find1), unique(find2)))
length(intersect(unique(find1_anno), unique(find2_anno)))

unique(find1)[order(unique(find1))]
unique(find2)[order(unique(find2))]

unique(find1_anno)[order(unique(find1_anno))]
unique(find2_anno)[order(unique(find2_anno))]








thresholds <- all_seed/100

thresholds <- 1:10

library(ggplot2)
library(dplyr)
library(tidyr)


ghost_counts <- as.numeric(table(number_gk_eur))
anno_counts  <- as.numeric(table(number_annogk_eur))

# thresholds <- all_seed/100

ghost_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(ghost_counts >= t)),
  Method = "GhostKnockoff"
)

anno_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(anno_counts >= t)),
  Method = "AnnoGK"
)

plot_df <- rbind(ghost_summary, anno_summary)

p1_thresh <- ggplot(plot_df, aes(x = factor(threshold), y = count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "", #"Number of signals detected ≥ threshold (EUR GWAS)",
    x = "Detection frequency threshold",
    y = "Number of signals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.title = element_blank()
  )

p1_thresh




ghost_counts <- as.numeric(table(number_gk_eas))
anno_counts  <- as.numeric(table(number_annogk_eas))

# thresholds <- all_seed/100

ghost_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(ghost_counts >= t)),
  Method = "GhostKnockoff"
)

anno_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(anno_counts >= t)),
  Method = "AnnoGK"
)

plot_df <- rbind(ghost_summary, anno_summary)

p2_thresh <- ggplot(plot_df, aes(x = factor(threshold), y = count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "", #"Number of signals detected ≥ threshold (EAS GWAS)",
    x = "Detection frequency threshold",
    y = "Number of signals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.title = element_blank()
  )

p2_thresh


library(patchwork)

g = (p1_thresh + theme(legend.position = "none")) | p2_thresh

g












ghost_counts <- as.numeric(table(find1))
anno_counts  <- as.numeric(table(find1_anno))

# thresholds <- all_seed/100

ghost_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(ghost_counts >= t)),
  Method = "GhostKnockoff"
)

anno_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(anno_counts >= t)),
  Method = "AnnoGK"
)

plot_df <- rbind(ghost_summary, anno_summary)

p1_thresh <- ggplot(plot_df, aes(x = factor(threshold), y = count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "", #"Number of signals detected ≥ threshold (EUR GWAS)",
    x = "Detection frequency threshold",
    y = "Number of signals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.title = element_blank()
  )

p1_thresh




ghost_counts <- as.numeric(table(find2))
anno_counts  <- as.numeric(table(find2_anno))

# thresholds <- all_seed/100

ghost_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(ghost_counts >= t)),
  Method = "GhostKnockoff"
)

anno_summary <- data.frame(
  threshold = thresholds/max(thresholds),
  count = sapply(thresholds, function(t) sum(anno_counts >= t)),
  Method = "AnnoGK"
)

plot_df <- rbind(ghost_summary, anno_summary)

p2_thresh <- ggplot(plot_df, aes(x = factor(threshold), y = count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "", #"Number of signals detected ≥ threshold (EAS GWAS)",
    x = "Detection frequency threshold",
    y = "Number of signals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.title = element_blank()
  )

p2_thresh


library(patchwork)

g = (p1_thresh + theme(legend.position = "none")) | p2_thresh

g





########### compare 10 groups:



list_find1 <- list()
list_find1_anno <- list()
list_find2 <- list()
list_find2_anno <- list()

# number_gk_eur <- list()
# number_annogk_eur <- list()
# number_gk_eas <- list()
# number_annogk_eas <- list()

#

all_seed = 100*(1:100)
# c(1,12,10,10000,123456, 100, 1000, 123, 1234, 12345)
for(seed in all_seed){
  region_1 <- numeric()
  region_anno_1 <- numeric()
  region_2 <- numeric()
  region_anno_2 <- numeric()
  lambdas1 <- list()
  lambdas2 <- list()

  current_find1 <- character()
  current_find2 <- character()
  current_find1_anno <- character()
  current_find2_anno <- character()


  threshold = 0.1
  sum1 = 0
  sum2 = 0
  chrid = '1'
  M = 1
  # seed = '12345' for dss


  chrlist = (1:22)[c(-21)]
  for(chrid in chrlist){
    print(paste0("chrid", chrid))
    ancestry0 = 'EUR'
    ancestry1 = c('EAS',"AFR","LAT")
    # ancestry1 = c('EAS')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("scz_baseline_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)

    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)

    for(i in 1:n_region){
      # print(i)
      snp = result$snp[[i]]
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
        current_find1 <- append(current_find1, paste(chrid, i, sep = '_'))
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
        current_find1_anno <- append(current_find1_anno, paste(chrid, i, sep = '_'))
      }
    }


    ancestry0 = 'EAS'
    ancestry1 = c('EUR',"AFR","LAT")
    # ancestry1 = c('EUR')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("scz_baseline_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)
    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)


    for(i in 1:n_region){
      # print(i)
      snp = result$snp[[i]]
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
        current_find2 <- append(current_find2, paste(chrid, i, sep = '_'))
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
        current_find2_anno <- append(current_find2_anno, paste(chrid, i, sep = '_'))
      }
    }

  }

  nrow(region_1)
  nrow(region_2)

  nrow(region_anno_1)
  nrow(region_anno_2)

  list_find1[[as.character(seed)]] <- current_find1
  list_find2[[as.character(seed)]] <- current_find2
  list_find1_anno[[as.character(seed)]] <- current_find1_anno
  list_find2_anno[[as.character(seed)]] <- current_find2_anno
}






############# compare the reproduce rate and Jaccard score between different ancestries

all_seeds <- as.character(100 * (1:100))
n_blocks <- 10
seeds_per_block <- 10

results_df <- data.frame(
  Block = 1:n_blocks,
  RepRate_Baseline = numeric(n_blocks),
  RepRate_AnnoGK = numeric(n_blocks),
  Jaccard_Baseline = numeric(n_blocks),
  Jaccard_AnnoGK = numeric(n_blocks)
)

for (b in 1:n_blocks) {
  # obtain the result for the certain group
  start_idx <- (b - 1) * seeds_per_block + 1
  end_idx <- b * seeds_per_block
  block_seeds <- all_seeds[start_idx:end_idx]

  # Baseline
  eur_base <- unique(unlist(list_find1[block_seeds]))
  eas_base <- unique(unlist(list_find2[block_seeds]))

  # AnnoGK
  eur_anno <- unique(unlist(list_find1_anno[block_seeds]))
  eas_anno <- unique(unlist(list_find2_anno[block_seeds]))

  # length(A ∩ B) / length(A)
  # ---------------------------------------------------------
  results_df$RepRate_Baseline[b] <- length(intersect(eur_base, eas_base)) / length(eur_base)
  results_df$RepRate_AnnoGK[b]   <- length(intersect(eur_anno, eas_anno)) / length(eur_anno)

  # length(A ∩ B) / length(A ∪ B)
  # ---------------------------------------------------------
  results_df$Jaccard_Baseline[b] <- length(intersect(eur_base, eas_base)) / length(union(eur_base, eas_base))
  results_df$Jaccard_AnnoGK[b]   <- length(intersect(eur_anno, eas_anno)) / length(union(eur_anno, eas_anno))
}

print("10 Blocks result:")
print(results_df)


cat("\n--- Replication Rate (EUR -> EAS) ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(results_df$RepRate_Baseline), sd(results_df$RepRate_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(results_df$RepRate_AnnoGK), sd(results_df$RepRate_AnnoGK)))

cat("\n--- Jaccard Index ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(results_df$Jaccard_Baseline), sd(results_df$Jaccard_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(results_df$Jaccard_AnnoGK), sd(results_df$Jaccard_AnnoGK)))


######### compare the reproduce rate among groups:

### eur:

all_seeds <- as.character(100 * (1:100))
n_blocks <- 10
seeds_per_block <- 10
freq_threshold = 0

# 1.
blocks_eur_base <- list()
blocks_eur_anno <- list()

for (b in 1:n_blocks) {
  start_idx <- (b - 1) * seeds_per_block + 1
  end_idx <- b * seeds_per_block
  print(paste0(start_idx, ":", end_idx))
  block_seeds <- all_seeds[start_idx:end_idx]

  base_all_hits <- unlist(list_find1[block_seeds])
  base_freq <- table(base_all_hits)
  blocks_eur_base[[b]] <- names(base_freq[base_freq > freq_threshold])

  anno_all_hits <- unlist(list_find1_anno[block_seeds])
  anno_freq <- table(anno_all_hits)
  blocks_eur_anno[[b]] <- names(anno_freq[anno_freq > freq_threshold])
}

# 2.calculate for C(10, 2) = 45 pairs
n_pairs <- (n_blocks * (n_blocks - 1)) / 2
pair_results <- data.frame(
  Pair_ID = character(n_pairs),
  RepRate_Baseline = numeric(n_pairs),
  RepRate_AnnoGK = numeric(n_pairs),
  Jaccard_Baseline = numeric(n_pairs),
  Jaccard_AnnoGK = numeric(n_pairs),
  stringsAsFactors = FALSE
)

counter <- 1
for (i in 1:(n_blocks - 1)) {
  for (j in (i + 1):n_blocks) {

    set_i_base <- blocks_eur_base[[i]]
    set_j_base <- blocks_eur_base[[j]]

    set_i_anno <- blocks_eur_anno[[i]]
    set_j_anno <- blocks_eur_anno[[j]]

    inter_base <- length(intersect(set_i_base, set_j_base))
    pair_results$RepRate_Baseline[counter] <- inter_base / mean(c(length(set_i_base), length(set_j_base)))
    pair_results$Jaccard_Baseline[counter] <- inter_base / length(union(set_i_base, set_j_base))

    inter_anno <- length(intersect(set_i_anno, set_j_anno))
    pair_results$RepRate_AnnoGK[counter] <- inter_anno / mean(c(length(set_i_anno), length(set_j_anno)))
    pair_results$Jaccard_AnnoGK[counter] <- inter_anno / length(union(set_i_anno, set_j_anno))

    pair_results$Pair_ID[counter] <- paste0("Block_", i, "_vs_", j)
    counter <- counter + 1
  }
}

# 3. results
cat("\n=== Between-Block Stability ===\n")
cat(sprintf("Total pairs: %d\n", n_pairs))

cat("\n--- Replication Rate (Overlap / Mean Size) ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(pair_results$RepRate_Baseline), sd(pair_results$RepRate_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(pair_results$RepRate_AnnoGK), sd(pair_results$RepRate_AnnoGK)))

cat("\n--- Jaccard Index ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(pair_results$Jaccard_Baseline), sd(pair_results$Jaccard_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(pair_results$Jaccard_AnnoGK), sd(pair_results$Jaccard_AnnoGK)))



### eas:


all_seeds <- as.character(100 * (1:100))
n_blocks <- 10
seeds_per_block <- 10
freq_threshold = 0

# 1.
blocks_eas_base <- list()
blocks_eas_anno <- list()

for (b in 1:n_blocks) {
  start_idx <- (b - 1) * seeds_per_block + 1
  end_idx <- b * seeds_per_block
  print(paste0(start_idx, ":", end_idx))
  block_seeds <- all_seeds[start_idx:end_idx]

  base_all_hits <- unlist(list_find2[block_seeds])
  base_freq <- table(base_all_hits)
  blocks_eas_base[[b]] <- names(base_freq[base_freq > freq_threshold])

  anno_all_hits <- unlist(list_find2_anno[block_seeds])
  anno_freq <- table(anno_all_hits)
  blocks_eas_anno[[b]] <- names(anno_freq[anno_freq > freq_threshold])
}

# 2.
n_pairs <- (n_blocks * (n_blocks - 1)) / 2
pair_results <- data.frame(
  Pair_ID = character(n_pairs),
  RepRate_Baseline = numeric(n_pairs),
  RepRate_AnnoGK = numeric(n_pairs),
  Jaccard_Baseline = numeric(n_pairs),
  Jaccard_AnnoGK = numeric(n_pairs),
  stringsAsFactors = FALSE
)

counter <- 1
for (i in 1:(n_blocks - 1)) {
  for (j in (i + 1):n_blocks) {

    set_i_base <- blocks_eas_base[[i]]
    set_j_base <- blocks_eas_base[[j]]

    set_i_anno <- blocks_eas_anno[[i]]
    set_j_anno <- blocks_eas_anno[[j]]

    inter_base <- length(intersect(set_i_base, set_j_base))
    pair_results$RepRate_Baseline[counter] <- inter_base / mean(c(length(set_i_base), length(set_j_base)))
    pair_results$Jaccard_Baseline[counter] <- inter_base / length(union(set_i_base, set_j_base))

    inter_anno <- length(intersect(set_i_anno, set_j_anno))
    pair_results$RepRate_AnnoGK[counter] <- inter_anno / mean(c(length(set_i_anno), length(set_j_anno)))
    pair_results$Jaccard_AnnoGK[counter] <- inter_anno / length(union(set_i_anno, set_j_anno))

    pair_results$Pair_ID[counter] <- paste0("Block_", i, "_vs_", j)
    counter <- counter + 1
  }
}

# 3.
cat("\n=== Between-Block Stability ===\n")
cat(sprintf("Total pairs: %d\n", n_pairs))

cat("\n--- Replication Rate (Overlap / Mean Size) ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(pair_results$RepRate_Baseline), sd(pair_results$RepRate_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(pair_results$RepRate_AnnoGK), sd(pair_results$RepRate_AnnoGK)))

cat("\n--- Jaccard Index ---\n")
cat(sprintf("Baseline: %.3f ± %.3f\n", mean(pair_results$Jaccard_Baseline), sd(pair_results$Jaccard_Baseline)))
cat(sprintf("AnnoGK:   %.3f ± %.3f\n", mean(pair_results$Jaccard_AnnoGK), sd(pair_results$Jaccard_AnnoGK)))
