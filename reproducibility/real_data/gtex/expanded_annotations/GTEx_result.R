# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)


source(repro_path("annokn_source", "KF_anno.R"))
source(repro_path("annogk_source", "GK_anno.R"))
source(repro_path("annogk_source", "GhostKnockoff.R"))
packageVersion("Matrix") 


sum_all_ct <- data.frame(matrix(nrow = 22, ncol = 0))
number_gene_ct <- data.frame(matrix(nrow = 22, ncol = 0))
find1_all <- list()
find2_all <- list()
find3_all <- list() # adakn



allcelltype = c('Muscle_Skeletal', 'Whole_Blood',"Skin_Sun_Exposed_Lower_leg","Artery_Tibial",
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")

#, "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic",
# "Lung", "Esophagus_Mucosa"

type = 1
alpha = 0.1
pc = 5
lnc = FALSE



library(data.table)
annot_all = fread(paste0(repro_path("ldsc_annotations", "baselineLD."), 1, ".annot.gz"))

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

a_dim <- length(cats) - 1




lambdas_all <- list()


for(celltype in allcelltype){
  
  lambda_s <- numeric()
  
  number_gene <- numeric()
  
  for(chrid in as.character(1:22)){
    print(celltype)
    print(chrid)
    library(data.table)
    library(genio)
    library(tidyr)
    library(dplyr)
    library(bigstatsr)
    
    
    bedNA <- function(bed1){
      for(j in 1:ncol(bed1)){
        temp <- bed1[,j]
        temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
        bed1[,j] <- temp
        #print(j)
      }
      return(bed1)
    }
    
    
    
    gexp_path = paste0(repro_path("gtex_expression", "adjusted_expr_age_sex_"),pc,'genopcs_5peers')
    
    directory_path <- paste(gexp_path, "/chr",chrid, sep="")
    # List all folders in the specified directory
    folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
    genes <- sub("\\..*", "", folder_names)
    
    # library(biomaRt)
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
        next  # Skip this iteration and move to the next one
      }
      remain <- append(remain, j)
    }
    
    risk_gene <- risk_gene[remain,]
    print(dim(risk_gene))
    number_gene <- append(number_gene, nrow(risk_gene))
  }
  
  
  
  sum_all <- numeric()
  sum1_all = c()
  sum2_all = c()
  sum3_all = c()
  find1 <- c()
  find2 <- c()
  find3 <- c()
  for(chrid in 1:22){
    # print(chrid)
    sum1 = c()
    sum2 = c()
    sum3 = c()
    path = paste0(repro_path("gtex_annotation_results", "chr_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
    
    if (!file.exists(path)) {
      message("not exist: ", path)
      break
    }
    load(path)
    
    # print(length(all_result))
    
    sum_all <- append(sum_all, length(all_result))
    
    for(gene_i in 1:length(all_result)){
      # print(gene_i)
      adakn_result = all_result[[gene_i]]
      
      W1 = adakn_result$ghostknockoff
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      # print("rej1:")
      # print(rej1)
      if(length(rej1 > 0)){
        sum1 = append(sum1, gene_i)
        find1 <- append(find1, paste(chrid, gene_i))
      }
      result_raw <- adakn_result$knockoff_anno
      
      
      if(length(unlist(result_raw$lambda_s)) < a_dim){
        print(paste("less than", a_dim,"annotation at celltype", celltype, "chr", adakn_result$chr, "gene", adakn_result$gene_i))
        # print(celltype)
        # print(adakn_result$chr)
        # print(adakn_result$gene_i)
        b = apply(adakn_result$R, 2, var)
        a = numeric(a_dim)
        a[b!=0] = unlist(result_raw$lambda_s)
        a[b==0] = 0
        lambda_s <- rbind(lambda_s, a)
      }else{
        lambda_s <- rbind(lambda_s, unlist(result_raw$lambda_s))
      }
      
      
      beta = result_raw$beta
      p = length(beta)/2
      T = abs(beta[1:p])
      T_tilde = abs(beta[(p+1):(2*p)])
      T_max = pmax(T,T_tilde)
      W = T - T_tilde
      # hist(W, breaks = 50)
      tau = knockoff.threshold(W,fdr = alpha,offset = 1)
      rej2 = as.numeric(which(W >= tau))
      
      # print("rej3:")
      # print(rej3)
      if(length(rej2 > 0)){
        sum2 = append(sum2, gene_i)
        find2 <- append(find2, paste(adakn_result$chr, adakn_result$gene_i))
      }
      
      result_raw = adakn_result$adakn
      
      if (!is.null(result_raw)){
        if(alpha == 0.2){
          rej3 <- result_raw$rejs[[1]]
        }else if(alpha == 0.1){
          rej3 <- result_raw$rejs[[2]]
        }else if(alpha == 0.05){
          rej3 <- result_raw$rejs[[3]]
        }
        if(length(rej3 > 0)){
          sum3 = append(sum3, gene_i)
          find3 <- append(find3, paste(adakn_result$chr, adakn_result$gene_i))
        }
      } 
      
      
    }
    sum1_all <- append(sum1_all, length(sum1))
    sum2_all <- append(sum2_all, length(sum2))
    sum3_all <- append(sum3_all, length(sum3))
  }
  
  cbind(number_gene, sum_all)
  
  print(length(find1))
  print(length(find2))
  
  sum_all_ct[[celltype]] = sum_all
  number_gene_ct[[celltype]] = number_gene
  
  find1_all <- append(find1_all, list(find1))
  find2_all <- append(find2_all, list(find2))
  find3_all <- append(find3_all, list(find3))
  
  lambdas_all <- append(lambdas_all, list(lambda_s))
  
}




sum_all_ct_0 = sum_all_ct
number_gene_ct_0 = number_gene_ct

number_gene_ct - sum_all_ct
# a = numeric()
# for(adakn_result in all_result){
#   a = append(a, adakn_result$gene_i)
# }


improve_prop_0 = (sapply(find2_all, length) - sapply(find1_all, length))/sapply(find1_all, length)
