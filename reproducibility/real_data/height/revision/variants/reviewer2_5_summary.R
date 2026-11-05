if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)



########## check the risk region for different CHRs

s = numeric()

for(chrid in 1:22){
  path = paste0(repro_path("height_single_regions", ""),"EUR","_chr_", chrid,".txt")
  risk_region <- read.table(path)
  if(any(is.na(risk_region))){
    s = append(s, 0)
    next
  }
  print(dim(risk_region))
  s = append(s, nrow(risk_region))
}


chr <- 1:length(s)
df <- data.frame(Chromosome = factor(chr), Count = s)

library(ggplot2)

g <- ggplot(df, aes(x = Chromosome, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Chromosome",
    y = "Region Count"
  )

g
# ggsave(repro_path("plots", "annoGk/Height_risk_region.pdf"),
#        g,
#        width = 7, height = 4, units = "in",
#        bg = "white", device = cairo_pdf)

########################

region_1 <- numeric()
region_anno_1 <- numeric()
region_2 <- numeric()
region_anno_2 <- numeric()
lambdas1 <- numeric()
lambdas2 <- numeric()

threshold = 0.1
sum1 = 0
sum2 = 0
chrid = '5'
M = 1
seed = '1000'
# seed = '12345' for dss

number_gk_eur <- 0
number_annogk_eur <- 0

number_gk_afr <- 0
number_annogk_afr <- 0


chrlist = (1:22)
for(chrid in chrlist){
  print(paste0("chrid", chrid))
  ancestry0 = 'EUR'
  ancestry1 = c('AFR')
  other = paste(ancestry1, collapse = "_")
  path = paste0(repro_path("height_single_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)

  print(length(result$risk_region))

  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)



  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      number_gk_eur <- number_gk_eur + sum(a <= threshold)
      region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      number_annogk_eur <- number_annogk_eur + sum(b <= threshold)
      region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas1 <- append(lambdas1, result$lambda_s[[i]])
    }

  }


  ancestry0 = 'AFR'
  ancestry1 = c('EUR')
  other = paste(ancestry1, collapse = "_")
  path = paste0(repro_path("height_single_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  print(length(result$risk_region))

  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)


  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      number_gk_afr <- number_gk_afr + sum(a <= threshold)
      region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      number_annogk_afr <- number_annogk_afr + sum(b <= threshold)
      region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas2 <- append(lambdas2, result$lambda_s[[i]])
    }

  }

}

nrow(region_1)
nrow(region_2)

nrow(region_anno_1)
nrow(region_anno_2)


df1 <- as.data.frame(region_1)
df2 <- as.data.frame(region_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common


df1 <- as.data.frame(region_anno_1)
df2 <- as.data.frame(region_anno_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common



########## Jaccard index

# 1.GhostKnockoff

df_gk_1 <- as.data.frame(region_1)
df_gk_2 <- as.data.frame(region_2)
common_gk <- merge(df_gk_1, df_gk_2)

intersect_gk <- nrow(common_gk)
union_gk <- nrow(df_gk_1) + nrow(df_gk_2) - intersect_gk
jaccard_gk <- intersect_gk / union_gk

# 2.AnnoGK
df_anno_1 <- as.data.frame(region_anno_1)
df_anno_2 <- as.data.frame(region_anno_2)
common_anno <- merge(df_anno_1, df_anno_2)

intersect_anno <- nrow(common_anno)
union_anno <- nrow(df_anno_1) + nrow(df_anno_2) - intersect_anno
jaccard_anno <- intersect_anno / union_anno

# 3. result
cat("========================================\n")
cat("GhostKnockoff Jaccard Index :", round(jaccard_gk, 4), "\n")
cat("AnnoGK Jaccard Index        :", round(jaccard_anno, 4), "\n")
cat("========================================\n")








library(ggplot2)
df <- data.frame(lambdas = lambdas1)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,
    fill = "#0072B2",
    color = "white",
    alpha = 0.8
  ) +
  labs(
    title = "",#expression("Distribution of " * lambda[l] * " of AFR p-values for EUR GWAS"),
    x = expression(lambda[l]),
    y = "Frequency (Count)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot

# ggsave(repro_path("plots", "annoGk/Height_lambda.pdf"),
#        histogram_plot,
#        width = 8, height = 5, units = "in",
#        bg = "white", device = cairo_pdf)


df <- data.frame(lambdas = lambdas2)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,
    fill = "#0072B2",
    color = "white",
    alpha = 0.8
  ) +
  labs(
    title = expression("Distribution of " * lambda[l] * " of EUR p-values for AFR GWAS"),
    x = expression(lambda[l]),
    y = "Frequency (Count)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot




library(VennDiagram)
library(grid)

df1 <- as.data.frame(region_1)
df2 <- as.data.frame(region_anno_1)
set1 <- apply(df1, 1, paste, collapse = "_")
set2 <- apply(df2, 1, paste, collapse = "_")

n1 <- length(set1)
n2 <- length(set2)
n_common <- length(intersect(set1, set2))


# pdf(repro_path("plots", "annoGk/Height_eur_venn.pdf"),
#     width = 5.5, height = 5)


venn.plot <- draw.pairwise.venn(
  area1 = n1,
  area2 = n2,
  cross.area = n_common,
  category = c("GhostKnockoff", "AnnoGK"),
  fill = c("steelblue", "orange"),

  alpha = c(0.6, 0.6),
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

# grid.newpage()
# grid.draw(venn.plot)

# grid.text(
#   "Risk regions identified by GhostKnockoff and AnnoGK in European",
#  x = 0.5, y = 0.98,
#  gp = gpar(fontsize = 16, fontface = "bold")
# )

# dev.off()




################# For AFR


library(VennDiagram)
library(grid)

df1 <- as.data.frame(region_2)
df2 <- as.data.frame(region_anno_2)
set1 <- apply(df1, 1, paste, collapse = "_")
set2 <- apply(df2, 1, paste, collapse = "_")

n1 <- length(set1)
n2 <- length(set2)
n_common <- length(intersect(set1, set2))


# pdf(repro_path("plots", "annoGk/Height_afr_venn.pdf"),
#     width = 5.5, height = 5)


venn.plot <- draw.pairwise.venn(
  area1 = n1,
  area2 = n2,
  cross.area = n_common,
  category = c("GhostKnockoff", "AnnoGK"),
  fill = c("steelblue", "orange"),

  alpha = c(0.6, 0.6),
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

# grid.newpage()
# grid.draw(venn.plot)

# grid.text(
#   "Risk regions identified by GhostKnockoff and AnnoGK in European",
#  x = 0.5, y = 0.98,
#  gp = gpar(fontsize = 16, fontface = "bold")
# )

# dev.off()



#################




region_1 <- numeric()
region_anno_1 <- numeric()
region_2 <- numeric()
region_anno_2 <- numeric()
lambdas1 <- numeric()
lambdas2 <- numeric()

snp_1 <- numeric()
snp_anno_1 <- numeric()
snp_2 <- numeric()
snp_anno_2 <- numeric()

threshold = 0.1
sum1 = 0
sum2 = 0
chrid = '5'
M = 1
seed = '3000'
# seed = '12345' for dss
all_seed = 100 * (1:10)

for(seed in all_seed){
  number_gk_eur <- 0
  number_annogk_eur <- 0

  number_gk_afr <- 0
  number_annogk_afr <- 0


  chrlist = (1:22)
  for(chrid in chrlist){
    print(paste0("chrid", chrid))
    ancestry0 = 'EUR'
    ancestry1 = c('AFR')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("height_ten_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)

    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)

    snp = result$snp



    for(i in 1:n_region){
      # print(i)
      a = q_gk[[i]]
      snp_temp = snp[[i]]
      if(sum(a <= threshold) > 0){
        number_gk_eur <- number_gk_eur + sum(a <= threshold)
        region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
        snp_1 <- append(snp_1, snp_temp$id[a <= threshold])
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        number_annogk_eur <- number_annogk_eur + sum(b <= threshold)
        region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas1 <- append(lambdas1, result$lambda_s[[i]])
        snp_anno_1 <- append(snp_anno_1, snp_temp$id[b <= threshold])
      }

    }


    ancestry0 = 'AFR'
    ancestry1 = c('EUR')
    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("height_ten_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)
    print(length(result$risk_region))

    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    snp = result$snp
    n_region = length(q_gk)


    for(i in 1:n_region){
      # print(i)
      snp_temp = snp[[i]]
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        number_gk_afr <- number_gk_afr + sum(a <= threshold)
        region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
        snp_2 <- append(snp_2, snp_temp$id[a <= threshold])
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        number_annogk_afr <- number_annogk_afr + sum(b <= threshold)
        region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas2 <- append(lambdas2, result$lambda_s[[i]])
        snp_anno_2 <- append(snp_anno_2, snp_temp$id[b <= threshold])
      }

    }

  }
}


region_1 <- apply(region_1, 1, paste, collapse = "_")
region_2 <- apply(region_2, 1, paste, collapse = "_")

region_anno_1 <- apply(region_anno_1, 1, paste, collapse = "_")
region_anno_2 <- apply(region_anno_2, 1, paste, collapse = "_")

length(unique(region_1))
length(unique(region_2))
length(unique(region_anno_1))
length(unique(region_anno_2))


common_rows <- intersect(unique(region_1), unique(region_2))
n_common <- length(common_rows)
n_common


common_rows <- intersect(unique(region_anno_1), unique(region_anno_2))
n_common <- length(common_rows)
n_common


length(unique(snp_1))
length(unique(snp_2))
length(unique(snp_anno_1))
length(unique(snp_anno_2))






###### jaccard index

# 1.GhostKnockoff

df_gk_1 <- as.data.frame(unique(region_1))
df_gk_2 <- as.data.frame(unique(region_2))
common_gk <- merge(df_gk_1, df_gk_2)

intersect_gk <- nrow(common_gk)
union_gk <- nrow(df_gk_1) + nrow(df_gk_2) - intersect_gk
jaccard_gk <- intersect_gk / union_gk

# 2.AnnoGK
df_anno_1 <- as.data.frame(unique(region_anno_1))
df_anno_2 <- as.data.frame(unique(region_anno_2))
common_anno <- merge(df_anno_1, df_anno_2)

intersect_anno <- nrow(common_anno)
union_anno <- nrow(df_anno_1) + nrow(df_anno_2) - intersect_anno
jaccard_anno <- intersect_anno / union_anno

# 3. result
cat("========================================\n")
cat("GhostKnockoff Jaccard Index :", round(jaccard_gk, 4), "\n")
cat("AnnoGK Jaccard Index        :", round(jaccard_anno, 4), "\n")
cat("========================================\n")
