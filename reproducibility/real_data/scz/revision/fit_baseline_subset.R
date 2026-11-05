if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)




args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)

chrid = as.numeric(args[1])

# source(repro_path("annokn_source", "KF_anno.R"))
# source(repro_path("annogk_source", "GK_anno.R"))
# source(repro_path("annogk_source", "GhostKnockoff.R"))
packageVersion("Matrix")


library(AnnoGKR)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  remotes::install_version("tidyr", version = "1.3.1", upgrade = "never")
}
library(tidyr)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  remotes::install_version("dplyr", version = "1.1.4", upgrade = "never")
}
library(dplyr)
if (!requireNamespace("genio", quietly = TRUE)) {
  remotes::install_version("genio", version = "1.1.2", upgrade = "never")
}
library(genio)
if (!requireNamespace("bigstatsr", quietly = TRUE)) {
  remotes::install_version("bigstatsr", version = "1.6.2", upgrade = "never")
}
library(bigstatsr)
if (!requireNamespace("dendextend", quietly = TRUE)) {
  remotes::install_version("dendextend", version = "1.19.1", upgrade = "never")
}
library(dendextend)
if (!requireNamespace("susieR", quietly = TRUE)) {
  remotes::install_version("susieR", version = "0.14.2", upgrade = "never")
}
library(susieR)


bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


################ load the summary statistics
ancestry0 = 'EAS'
ancestry1 = c('EUR',"AFR","LAT")


ancestry = ancestry0
path = paste0(repro_path("scz_gwas", ""), ancestry, "_summ.tsv")
sum_s <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# we also checked that all SNPs are in (A,T,G,C)
colnames(sum_s) <- c("Chromsome","SNPID", "Position", "EffectAllele", "NonEffectAllele", "NonEffectAF_case",
                     "NonEffectAF_control", "INFO", "Beta", "SE","Pval","N_case","N_control","Neff")
sum_s[[ancestry]] <- sum_s$Pval
print(length(unique(sum_s$SNPID)) == nrow(sum_s))

# ,'HIS','SAS','EAS','AFA'
for(ancestry in ancestry1){
  path = paste0(repro_path("scz_gwas", ""), ancestry, "_summ.tsv")
  temp <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  colnames(temp) <- c("Chromsome","SNPID", "Position", "EffectAllele", "NonEffectAllele", "NonEffectAF_case",
                      "NonEffectAF_control", "INFO", "Beta", "SE","Pval","N_case","N_control","Neff")

  message(ancestry, ": unique SNP? ", length(unique(temp$SNPID)) == nrow(temp))
  index <- match(sum_s$SNP, temp$SNP)
  print(table(is.na(index)))
  print(table(temp$SNP[index] == sum_s$SNP))
  message(ancestry, ": matched = ", sum(!is.na(index)), "/", length(index))
  sum_s[[ancestry]] <- temp$Pval[index]
}




###############


M = 1
seed = 12345

# 100, 1000, 123, 1234, 12345

all_seed = 100*(51:100)

for(M in c(1,3)){
  for(seed in all_seed){
    sums_chr = sum_s[sum_s$Chromsome == chrid,]
    sums_chr <- sums_chr[order(sums_chr$Position),]

    ################ load the reference panel from 1KG

    # ancestry = 'EUR'
    if(ancestry0 == 'EUR'){
      path <- paste(repro_path("cross_ancestry_plink", "bychr_"), ancestry0, "/1KG_chr", chrid,sep='')
    }else if(ancestry0 == 'EAS'){
      path <- paste(repro_path("cross_ancestry_plink", "bychr_"), ancestry0, "/1KG_chr", chrid,sep='')
    }
    ref_chr <- read_plink(path)

    ############### load the annotation data

    # library(data.table)
    # cats = c("SNP","Coding_UCSC", "TSS_Hoffman", "Promoter_UCSC", "Intron_UCSC", "UTR_3_UCSC", "UTR_5_UCSC")
    # annot = fread(paste0(repro_path("archived_ldsc_annotations", "baselineLD."), chrid, ".annot.gz"))
    # annot = annot[, ..cats]


    ############### find the common part of three datasets:

    common_id <- intersect(ref_chr$bim$pos, sums_chr$Position)
    index1 = match(common_id, ref_chr$bim$pos)
    index2 = match(common_id, sums_chr$Position)

    bim = ref_chr$bim[index1, ]
    X = ref_chr$X[index1, ]
    table(rownames(X) == bim$id)
    sums_chr = sums_chr[index2, ]
    table(bim$pos == sums_chr$Position)


    ################### define risk regions

    # len_half = 250000
    # # risk_id <- which(sums_chr$Pval < 5 * 10^-8)
    # risk_id <- which(sums_chr$EUR < 5 * 10^-8)
    # pos_id = unlist(sums_chr$Position[risk_id])
    # pos_diff <- diff(pos_id)
    # where_sep = as.numeric(c(0, which(pos_diff > 2 * len_half)))
    # risk_region = numeric()
    #
    #
    # if(length(where_sep) > 1){
    #   for(i in 1:(length(where_sep)-1)){
    #     risk_region <- rbind(risk_region, c(pos_id[where_sep[i] + 1] - len_half, pos_id[where_sep[i+1]] + len_half))
    #   }
    #   risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
    # }else if(length(where_sep) == 1){
    #   risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
    # }
    #
    # path = paste0(repro_path("archived_scz_regions", ""),ancestry0,"_chr_", chrid,".txt")
    # write.table(risk_region, path, row.names = FALSE, col.names = FALSE, quot = FALSE)
    #
    #### risk_region <- risk_region[2:nrow(risk_region),] remove the first region (MHC on chromosome 6)


    path = paste0(repro_path("scz_candidate_regions", ""),"EUR","_chr_", chrid,".txt")
    risk_region <- read.table(path)

    ################### address allele flipping

    q_gk <- list()
    q_gkanno <- list()
    lambda_s <- list()
    risk_region_final <- list()
    snp_information <- list()

    i = 1

    for(i in 1:nrow(risk_region)){
      print(paste("region:", i))

      set.seed(seed)

      start = risk_region[i, 1]
      end = risk_region[i, 2]


      region_id = which(bim$pos > start & bim$pos < end)
      X_mid = bedNA(t(X[region_id, ]))
      na_cols <- apply(X_mid, 2,function(x) all(is.na(x)))
      which(na_cols)
      same_cols <- apply(X_mid, 2,function(x) length(unique(x)))
      table(same_cols)

      region_id = region_id[which(same_cols > 1)]


      trait_region <- sums_chr[region_id, ]
      bim_region <- bim[region_id, ]
      X_region <- bedNA(t(X[region_id, ]))
      # annot_region <- annot[region_id, ]
      # table(colnames(X_region) == annot_region$SNP)
      # table(bim_region$id == annot_region$SNP)
      # table(trait_region$Position == bim_region$pos)


      z_scores = trait_region$Beta/trait_region$SE
      R = cor(X_region)
      nsample <- median(trait_region$Neff)
      Neff <- trait_region$Neff


      library(susieR)
      lambda = estimate_s_rss(z_scores, R, n = nsample)
      print(lambda)


      inv <- which(trait_region$EffectAllele==bim_region$ref & trait_region$NonEffectAllele==bim_region$alt)
      print(length(inv))

      trait_region$Beta[inv] <- -trait_region$Beta[inv]
      trait_region$EAF[inv] <- 1-trait_region$EAF[inv]
      a <- trait_region$EffectAllele[inv]
      trait_region$EffectAllele[inv] <- trait_region$NonEffectAllele[inv]
      trait_region$NonEffectAllele[inv] <- a

      print(table(trait_region$EffectAllele==bim_region$alt & trait_region$NonEffectAllele==bim_region$ref))


      falseid <- which((trait_region$EffectAllele==bim_region$alt & trait_region$NonEffectAllele==bim_region$ref)==F)

      print(paste("number of falseid:",length(falseid)))

      if(length(falseid)>0){
        print("removing false SNPs:")
        print(falseid)
        trait_region <- trait_region[-falseid,]
        bim_region <- bim_region[-falseid,]
        X_region <- X_region[,-falseid]
      }

      print(table(trait_region$EffectAllele==bim_region$alt & trait_region$NonEffectAllele==bim_region$ref))


      z_scores = trait_region$Beta/trait_region$SE

      R = cor(X_region)
      # library(susieR)
      lambda = estimate_s_rss(z_scores, R, n = nsample)
      print(lambda)



      ################### perform variable selection


      #### use hierarchical clustering to cluster SNPs

      # d <- as.dist(1 - R)  # higher distance = lower correlation
      R2 <- R^2
      d <- as.dist(1 - R2) # cluster SNPs with r² > 0.25
      hc <- hclust(d, method = "average")
      hc <- as.dendrogram(hc)
      clusters <- cutree(hc, h = 0.75)  # e.g., cluster SNPs with r² > 0.25
      length(table(clusters))
      mean(table(clusters)) # average size for each cluster


      ######### select representative SNP

      m_annot <- length(ancestry1)

      table(names(clusters) == bim_region$id)

      df <- data.frame(
        SNP = names(clusters),
        pvalues = trait_region$Pval,
        cluster = clusters,
        index = 1:nrow(trait_region)
      )

      representatives <- df %>%
        group_by(cluster) %>%
        slice_min(order_by = pvalues, n = 1, with_ties = FALSE)

      head(representatives)
      dim(representatives)
      index <- representatives$index[order(representatives$index)]
      X_final = X_region[, index]
      bim_final <- bim_region[index,]
      z_scores_final  = z_scores[index]
      # annot_final = data.frame(position = trait_region[index, 2])
      annot_final = trait_region[index, c(2,(ncol(trait_region) - m_annot + 1):ncol(trait_region))]
      # colnames(annot_final)[1] = 'Position'
      nsample= median(trait_region$Neff[index])
      Neff = trait_region$Neff[index]

      ##########
      X_final = scale(X_final)
      table(colnames(X_final) == bim_final$id)
      table(bim_final$id == annot_final$SNPID)
      LD = cor(X_final)
      p = nrow(LD)


      ############################# ghostknockoff

      fit.prelim <- GhostKnockoff.prelim(
        cor.G   = LD,
        M       = M,
        method  = "sdp"
      )

      GK1_lasso <- GhostKnockoff.fit(z_scores_final, Neff, fit.prelim, method='lasso')
      GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
      threshold = 0.1
      which(GK.filter$q <= threshold)
      a = GK.filter$q

      q_gk <- append(q_gk, list(a))


      ############################# GK-anno

      annot_final = annot_final[, 2:ncol(annot_final)]
      annot_final_df <- as.data.frame(annot_final)

      annot_final_df <- apply(annot_final_df, 2, function(col){
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        -log(col)
      })
      annot_final_df <- annot_final_df[, which(apply(annot_final_df, 2, var) > 0)]


      R_anno = scale(annot_final_df)
      GK1ps_anno = GK_anno_M(z_scores_final, R_anno, M, LD, nsample)
      # GK1ps_anno = GK_anno_dss(z_scores_final, R_anno, M, LD, Neff)

      # beta <- GK1ps_anno$beta_final
      # T_0<-abs(beta[1:p])
      # T_k<-abs(matrix(beta[-(1:p)],p,M))
      T_0 <- GK1ps_anno$T_0
      T_k <- GK1ps_anno$T_k

      GK.filter<-GhostKnockoff.filter(T_0,T_k)
      threshold = 0.1
      which(GK.filter$q <= threshold)
      b = GK.filter$q

      q_gkanno <- append(q_gkanno, list(b))

      lambda_s <- append(lambda_s, list(GK1ps_anno$lambda_s))

      risk_region_final <- append(risk_region_final, list(risk_region[i,]))
      snp_information <- append(snp_information, list(bim_final))


    }

    result <- list(chr = chrid,
                   risk_region = risk_region_final,
                   snp = snp_information,
                   q_gk = q_gk,
                   q_gkanno = q_gkanno,
                   lambda_s = lambda_s)

    other = paste(ancestry1, collapse = "_")
    path = paste0(repro_path("scz_baseline_results", ""), ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    save(result, file = path)



  }
}
