if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)



#### note: for type 1 and type 2, we used the lams <- exp(seq(log(lambda)-1, log(lambda)+2, 0.1)), with more breaks and candidate lambda_0
### for type 3, the annotations dimensional is around 100, we try to use the lams <- exp(seq(log(lambda)-1, log(lambda)+2, 0.2)) to include less lambda_0?
### for type 3 now, we still use lams <- exp(seq(log(lambda)-1, log(lambda)+2, 0.1)), as only a few cases may fail


args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)

chrid = as.numeric(args[1])

#source(repro_path("annokn_source", "KF_anno.R"))
#source(repro_path("annogk_source", "GK_anno.R"))
#source(repro_path("annogk_source", "GhostKnockoff.R"))
packageVersion("Matrix")


library(data.table)
if (!requireNamespace("genio", quietly = TRUE)) {
  remotes::install_version(
    "genio",
    version = "1.1.2",
    upgrade = "never")
}
library(genio)
library(tidyr)
library(dplyr)
if (!requireNamespace("bigstatsr", quietly = TRUE)) {
  remotes::install_version(
    "bigstatsr",
    version = "1.6.2",
    upgrade = "never"
  )
}
library(bigstatsr)
library(knockoff)
library(glmnet)
library(AnnoKn)
# library(biomaRt)


bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


type = 1
pc = 5
# chrid = 1
lnc = FALSE
celltype = 'Skin_Sun_Exposed_Lower_leg'
# "Brain_Cortex" "Lung" "Pancreas"
print(celltype)

gexp_path = paste0(repro_path("gtex_expression", "adjusted_expr_age_sex_"),pc,'genopcs_5peers')

directory_path <- paste(gexp_path, "/chr",chrid, sep="")
# List all folders in the specified directory
folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
genes <- sub("\\..*", "", folder_names)


# #################################   avoid access online resource
#
# options(timeout = 300)
#
# while (TRUE) {
#   ensembl <- try({
#     # useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#     useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
#   })
#   if (!isa(ensembl, "try-error"))
#     break
#   else{
#     Sys.sleep(3)
#     print(ensembl)
#     print("try ensembl again!")
#   }
# }
#
# print("useMart")
#
# while (TRUE) {
#   z <- try({
#     getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "transcript_biotype"),
#           "ensembl_gene_id", genes, mart = ensembl)
#   })
#   if (!isa(z, "try-error"))
#     break
#   else{
#     Sys.sleep(3)
#     print(z)
#     print("try z again!")
#   }
#
# }
#
# index <- match(z$ensembl_gene_id, genes)
# z$folder_names <- folder_names[index]
# risk_gene <- z
#
# ####### save result
#
# path = paste0(repro_path("gtex_gene_metadata", "chr_"), chrid, "_ct_", celltype, ".txt")
# # write.table(risk_gene, path,
# #             row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(risk_gene, path,
#             row.names = FALSE,
#             col.names = TRUE,
#             quote = FALSE,
#             sep = "\t",
#             na = "NA")
#
# ###########

#################################

path = paste0(repro_path("gtex_gene_coordinates", "chr_"), chrid, "_ct_", celltype, ".txt")
risk_gene = fread(path)



if(lnc == FALSE){
  risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding"),]
}else{
  risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding", "lncRNA"),]
}

risk_gene <- risk_gene[!duplicated(risk_gene$ensembl_gene_id),]


remain <- numeric()
for(j in 1:nrow(risk_gene)){
  gene_info <- risk_gene[j,]
  path <- paste(gexp_path, "/chr",
                chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")
  if (!file.exists(path)) {
    cat("File does not exist for iteration", j, "\n")
    next  # Skip this iteration and move to the next one
  }
  remain <- append(remain, j)
}

risk_gene <- risk_gene[remain,]



########### read in the genotype data and match to snp

maf = 0.01

if(maf == 0.05){
  maf_indi = 'maf5'
  maf_dir = 'bychr_maf0.05'
}else if(maf == 0.01){
  maf_indi = 'maf1'
  maf_dir = 'bychr_maf0.01'
}

ancestry = 'EUR'
path <- paste(repro_path("gtex_plink", ""),ancestry,"/", maf_dir, "/GTEX_chr", chrid,sep='')
gtex_chr <- read_plink(path)

# snp_mart <- useEnsembl(
#   biomart = "ENSEMBL_MART_SNP",
#   dataset = "hsapiens_snp"
# )

file_path = paste0(repro_path("ldsc_annotations", "baselineLD."), chrid, ".annot.gz")
# annot_all = fread(file_path)
annot_all <- as.data.table(read.table(gzfile(file_path), header = TRUE))

cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "Coding_UCSC", "TFBS_ENCODE")
if(type == 1){
  cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "TFBS_ENCODE")
}else if(type == 2){
  cats <- c(
    "SNP",
    # QTL
    "GTEx_eQTL_MaxCPP",
    "BLUEPRINT_H3K27acQTL_MaxCPP",
    "BLUEPRINT_H3K4me1QTL_MaxCPP",
    "BLUEPRINT_DNA_methylation_MaxCPP",

    # Promoter / TSS
    "Promoter_UCSC",
    "TSS_Hoffman",
    "PromoterFlanking_Hoffman",
    "Human_Promoter_Villar",
    "Ancient_Sequence_Age_Human_Promoter",
    "Human_Promoter_Villar_ExAC",

    # Enhancer / TFBS
    "Enhancer_Andersson",
    "Enhancer_Hoffman",
    "SuperEnhancer_Hnisz",
    "WeakEnhancer_Hoffman",
    "TFBS_ENCODE",
    "Human_Enhancer_Villar",
    "Ancient_Sequence_Age_Human_Enhancer",

    # Histone marks
    "H3K27ac_Hnisz",
    "H3K27ac_PGC2",
    "H3K4me1_peaks_Trynka",
    "H3K4me1_Trynka",
    "H3K4me3_peaks_Trynka",
    "H3K4me3_Trynka",
    "H3K9ac_peaks_Trynka",
    "H3K9ac_Trynka"
  )
}else if(type == 3){
  cats <- colnames(annot_all)[c(-1, -2, -4, -5)]
}
# cats = c("SNP","Coding_UCSC", "TSS_Hoffman", "Enhancer_Hoffman",
#          "Promoter_UCSC", "H3K27ac_Hnisz", "H3K4me1_peaks_Trynka")
# annot = fread(paste0(repro_path("archived_ldsc_annotations", "baselineLD."), chrid, ".annot.gz"))
annot_all = annot_all[, ..cats]

######################################################################################
# path = paste0(repro_path("gtex_variant_coordinates", "chr"),chrid,".txt")
# ref_hg38 <- read.table(path, header = FALSE)
# colnames(ref_hg38) <- c("chr", "start","end","snp")
#
# ref_hg38_selected <- ref_hg38[ref_hg38$end %in% gtex_chr$bim$pos | ref_hg38$snp %in% annot_all$SNP,]
# path = paste0(repro_path("gtex_variant_coordinates", "chr"),chrid,"_selected.txt")
# write.table(ref_hg38_selected, file = path, row.names = FALSE,
#             col.names = TRUE, quote = FALSE)
######################################################################################

path = paste0(repro_path("gtex_variant_coordinates", "chr"),chrid,"_selected.txt")
ref_hg38 <- read.table(path, header = TRUE)


####################
# filtered_ref <- ref_hg38 %>%
#   group_by(end) %>%
#   slice(if (any(snp %in% annot_all$SNP)) {
#     which(snp %in% annot_all$SNP)[1]
#   } else {
#     sample(1:n(), 1)
#   }) %>%
#   ungroup()

##############

M = 250000

gene_i = 1



# path = paste0(repro_path("archived_gtex_results", "chr_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
all_result <- list()
# save(all_result, file = path)

for(gene_i in 1:nrow(risk_gene)){

  print(paste(gene_i, "out of", nrow(risk_gene)))
  set.seed(1000)

  gene_info <- risk_gene[gene_i,]
  path <- paste(gexp_path, "/chr",
                chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")

  start = gene_info$start_position - M
  end = gene_info$end_position + M

  # check
  if(j==1){
    print(paste("read gene expression data from", path))
  }

  expr <- fread(path)
  expr = as.data.frame(expr)
  short_ids <- sapply(strsplit(expr$individual, "-"), function(x) paste(x[1:2], collapse = "-"))
  expr$individual = short_ids
  colnames(expr) = c("V1","V2")



  ######### snps that within this gene:
  index <- which(gtex_chr$bim$pos > start & gtex_chr$bim$pos < end)
  # dim(gtex_chr$X)
  X_gtex_region <- gtex_chr$X[index, ]

  individuals <- colnames(X_gtex_region)
  common_ind <- intersect(expr$V1, individuals)
  gidx <- match(common_ind, expr$V1)
  midx <- match(common_ind, individuals)

  expr <- expr[gidx,]
  y = expr$V2

  if(var(y) == 0){
    cat("variance of gene expression is 0:",
        "  chrid:", chrid,
        "  gene:", gene_i, "\n",
        file = paste0(repro_path("gtex_revision_logs", "log_"), celltype, ".txt"),
        append = TRUE)
    next
  }


  ####### remove constant SNPs
  X_region <- X_gtex_region[,midx]
  maf <- apply(X_region, 1, function(x){mean(x, na.rm = TRUE)})
  bad_snp = as.numeric(which(apply(X_region, 1, var)==0))

  if(length(bad_snp) > 0){
    print(paste0("gene:", j, "  bad_snp:", bad_snp))
    # snpindex = snpindex[-bad_snp]
    X_region <- X_region[-bad_snp,]
  }

  ######## cluster and select representative SNPs

  X_region <- t(bedNA(X_region))
  X_region <- scale(X_region)

  table(rownames(X_region) == expr$V1)

  ########## calculate the pvalue

  X_fbm <- as_FBM(X_region)
  res <- big_univLinReg(X_fbm, y.train = y)
  pvals <- predict(res, log10 = FALSE)

  ################


  R = cor(X_region)
  R2 <- R^2
  d <- as.dist(1 - R2)
  # d <- as.dist(1 - R)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- tryCatch(
    cutree(hc, h = 0.75),
    error = function(e) {
      message("average linkage failed, switching to single linkage...")
      hc <<- hclust(d, method = "single")
      plot(hc, labels = FALSE, main = "SNP clustering by LD (single linkage)")
      cutree(hc, h = 0.75)
    }
  )
  # clusters <- cutree(hc, h = 0.75)  # e.g., cluster SNPs with r² > 0.5

  length(table(clusters))
  mean(table(clusters)) # average size for each cluster

  if(length(table(clusters)) < 11){
    cat("number of SNP groups <= 10:",
        "  chrid:", chrid,
        "  gene:", gene_i, "\n",
        file = paste0(repro_path("gtex_revision_logs", "log_"), celltype, ".txt"),
        append = TRUE)
    next
  }


  df <- data.frame(
    SNP = names(clusters),
    pval = pvals,
    cluster = clusters,
    index = 1:length(pvals)
  )

  representatives <- df %>%
    group_by(cluster) %>%
    slice_min(order_by = pval, n = 1, with_ties = FALSE)

  head(representatives)
  dim(representatives)


  index <- representatives$index[order(representatives$index)]

  X_repre = X_region[, index]
  y_repre = (y - mean(y))/sd(y)


  print(dim(X_repre))



  #################### convert to rsid

  name_snp <- colnames(X_repre)
  df <- data.frame(name_snp = name_snp)
  df_split <- df %>%
    separate(name_snp, into = c("chr", "pos", "ref", "alt", "build"), sep = "_") %>%
    mutate(chr = gsub("^chr", "", chr),
           pos = as.integer(pos))

  head(df_split)
  # position <- paste(df_split$chr, df_split$pos, df_split$pos, sep = ":")


  length(intersect(df_split$pos, ref_hg38$end))
  index <- match(df_split$pos, ref_hg38$end)
  snp_name <- ref_hg38$snp[index]

  ############################################################
  ref_hg38_temp <- ref_hg38[ref_hg38$end %in% df_split$pos,]
  ref_hg38_temp <- ref_hg38_temp %>%
    group_by(end) %>%
    slice(if (any(snp %in% annot_all$SNP)) {
      which(snp %in% annot_all$SNP)[1]
    } else {
      sample(1:n(), 1)
    }) %>%
    ungroup()

  index <- match(df_split$pos, ref_hg38_temp$end)
  snp_name <- ref_hg38_temp$snp[index]
  ############################################################


  # is snp not exists in ref_hg38, which is very unlikely, name this SNP the position on this chromosome
  snp_name[is.na(snp_name)] <- df_split$pos[is.na(snp_name)]
  colnames(X_repre) <- snp_name


  ############## knockoff

  # set.seed(111)
  Xk = create.second_order(X_repre, method = 'sdp')

  # standardize the matrix
  Xk = scale(Xk)

  # run original knockoff
  p = ncol(X_repre)
  mdl = cv.glmnet(cbind(X_repre, Xk), y_repre, alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[, mdl$lambda == cvlambda]
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W1 = T-T_tilde
  as.numeric(W1)
  hist(W1, breaks = 50)

  alpha = 0.2
  tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
  rej1 = as.numeric(which(W1>=tau))
  print(rej1)



  ############## knockoff with annotation
  # "Enhancer_Andersson", "GTEx_eQTL_MaxCPP"

  length(intersect(colnames(X_repre), annot_all$SNP))
  index <- match(colnames(X_repre), annot_all$SNP)
  annot <- annot_all[index, ]
  annot0 <- annot
  annot <- annot[, -1]

  apply(annot, 2, table)

  R <- apply(annot, 2, function(col){
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })


  ############################### AdaKn ################################################
  # alphalist <- c(0.05, 0.1, 0.2, 0.3)
  # res_result = filter_randomForest(W1,R,alpha =alphalist,offset=1)
  # rej3 = res_result$rejs[[2]]
  # rej3

  ######################################################################################


  result_raw <- AnnoKn(X = X_repre, Xk = Xk, y = y_repre, attempts = c(0), R = R)

  beta = result_raw$beta
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W = T - T_tilde

  hist(W, breaks = 50)


  tau = knockoff.threshold(W,fdr = alpha,offset = 1)
  rej2 = as.numeric(which(W >= tau))

  result_raw$lambda_s

  # hist(W[W!=0],  breaks = 20)
  # hist(W1[W1!=0],  breaks = 20)
  # hist(result_raw$Prior_final, xlab = 'weight', main = 'Histogram of weights')


  print(result_raw$lambda_s)
  print(rej1)
  print(rej2)


  final_result = list(chr = chrid,
                      gene_i = gene_i,
                      ghostknockoff = W1,
                      knockoff_anno = result_raw,
                      snp_name = snp_name,
                      annot = annot0,
                      R = R,
                      gene_info = gene_info)

  # final_result = list(
  #                     ghostknockoff = W1,
  #                     # knockoff_anno = result_raw,
  #                     adakn = res_result,
  #                     snp_name = snp_name,
  #                     annot = annot0,
  #                     R = R,
  #                     gene_info = gene_info)
  all_result <- append(all_result, list(final_result))

}

path = paste0(repro_path("gtex_revision_results", "chr_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
save(all_result, file = path)
