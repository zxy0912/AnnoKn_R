# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)

get_positive_proportion <- function(W, n_bins = 10) {
  
  W <- W[is.finite(W) & W != 0]
  
  absW <- abs(W)
  
  breaks <- quantile(
    absW,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )
  
  breaks <- unique(breaks)
  
  bin <- cut(
    absW,
    breaks = breaks,
    include.lowest = TRUE
  )
  
  result <- aggregate(
    W,
    by = list(bin = bin),
    FUN = function(x) mean(x > 0)
  )
  
  colnames(result)[2] <- "positive_proportion"
  
  result
}



chrid = '1'
celltype = 'Muscle_Skeletal'
type = '1'
simulation = 'TRUE'
# simulation = ''
# Adipose_Subcutaneous  Muscle_Skeletal Whole_Blood


path = paste0(repro_path("gtex_null_results", "knockoff_chr_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
if(simulation == 'TRUE'){
  path = paste0(repro_path("gtex_null_results", "knockoff_chr_simu_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
}
load(path)

path = paste0(repro_path("gtex_null_results", "annokn_chr_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
if(simulation == 'TRUE'){
  path = paste0(repro_path("gtex_null_results", "annokn_chr_simu_"), chrid, "_ct_", celltype, "_type_", type, ".RData")
}
load(path)

ls()

library(ggplot2)

df_hist <- data.frame(
  W = c(W_Knockoff[W_Knockoff != 0],
        W_Annokn[W_Annokn != 0]),
  Method = c(
    rep("Knockoff", sum(W_Knockoff != 0)),
    rep("AnnoKn", sum(W_Annokn != 0))
  )
)

g = ggplot(df_hist, aes(x = W, fill = Method)) +
  geom_histogram(
    bins = 50,
    alpha = 0.5,
    position = "identity"
  ) +
  labs(
    x = "W",
    y = "Count",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  )


ggsave(repro_path("revision_plots", "GTEx4.1/null_distribution.pdf"),
       g, 
       width = 6, height = 4, units = "in", 
       bg = "white", device = cairo_pdf)



get_positive_proportion(W_Knockoff[W_Knockoff != 0])
get_positive_proportion(W_Annokn[W_Annokn != 0])

binom.test(
  sum(W_Knockoff > 0),
  sum(W_Knockoff != 0),
  p = 0.5
)



###################
library(ggplot2)

df_knockoff <- get_positive_proportion(W_Knockoff[W_Knockoff != 0])
df_annokn <- get_positive_proportion(W_Annokn[W_Annokn != 0])

df_knockoff$Method <- "Knockoff"
df_annokn$Method <- "AnnoKn"

df_prop <- rbind(
  df_knockoff[, c("bin", "positive_proportion", "Method")],
  df_annokn[, c("bin", "positive_proportion", "Method")]
)

g = ggplot(df_prop, aes(x = bin, y = positive_proportion, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed"
  ) +
  labs(
    x = "|W| bin",
    y = "Positive proportion",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(repro_path("revision_plots", "GTEx4.1/null_proportion.pdf"),
       g, 
       width = 8, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)

