# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)

N.effect = 10000 # sample size
n1 = 1000 # sample size used for clustering
# p = 200 # number of variable
h2e = 0.02 # heritability 0.02, 0.05, 0.1, 0.2
alphae = 0.3 # proportion of causal SNPs. 0.05, 0.1, 0.2
# rho = 0.5 # maximum correlation between SNPs
M = 1 # number of knockoff copies
chrid = '1'
regu = 0.2
width = "2MB"
if(width == '2MB'){
  p1 = 100
}else{
  p1 = 200
}

power1 <- power2 <- power3 <- power4 <- power5 <- power6 <- power7 <- numeric()
fdr1 <- fdr2 <- fdr3 <- fdr4 <- fdr5 <- fdr6 <- fdr7<- numeric()
lambda <- numeric()

for (iter in 1:100) {
  
  path = paste0(repro_path("true_genotype_legacy_results", "result_"), width, "_insample_iter_", iter, "_n_", N.effect, "_h2e_", h2e, "_alphae_", alphae, "_p1_", p1, "_M_", M, "_chr_", chrid, ".RData")
  load(path)
  
  power1 <- rbind(power1, result$power1)
  power2 <- rbind(power2, result$power2)
  power3 <- rbind(power3, result$power3)
  power4 <- rbind(power4, result$power4)
  power5 <- rbind(power5, result$power5)
  power6 <- rbind(power6, result$power6)
  power7 <- rbind(power7, result$power7)
  
  fdr1 <- rbind(fdr1, result$fdr1)
  fdr2 <- rbind(fdr2, result$fdr2)
  fdr3 <- rbind(fdr3, result$fdr3)
  fdr4 <- rbind(fdr4, result$fdr4)
  fdr5 <- rbind(fdr5, result$fdr5)
  fdr6 <- rbind(fdr6, result$fdr6)
  fdr7 <- rbind(fdr7, result$fdr7)
  
  lambda = append(lambda, result$lambda)
}

colMeans(power1, na.rm = TRUE)
colMeans(power2, na.rm = TRUE)
colMeans(power3, na.rm = TRUE)

colMeans(fdr1, na.rm = TRUE)
colMeans(fdr2, na.rm = TRUE)
colMeans(fdr3, na.rm = TRUE)

colMeans(power4, na.rm = TRUE)
colMeans(power5, na.rm = TRUE)
colMeans(power6, na.rm = TRUE)


colMeans(fdr4, na.rm = TRUE)
colMeans(fdr5, na.rm = TRUE)
colMeans(fdr6, na.rm = TRUE)








library(ggplot2)

methods <- c("Knockoff", "AnnoKn-simple", "AnnoKn", 
             "GhostKnockoff", "AnnoGK-simple", "AnnoGK", "AdaKn (RF)")


colors <- c("Knockoff" = "#E69F00",
            "AnnoKn" = "#FF0000",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#e7298a",
            "AnnoKn-simple" = "#66a61e",
            "AnnoGK-simple" = "#e6ab02",
            "AdaKn (RF)" = "#1b9e77")

shapes <- c("Knockoff" = 16,
            "AnnoKn" = 8,
            "GhostKnockoff" = 15,
            "AnnoGK" = 3,
            "AnnoKn-simple" = 5,
            "AnnoGK-simple" = 4,
            "AdaKn (RF)" = 17)

fdr_seq <- seq(0.4, 0.1, by = -0.05)
len = length(fdr_seq)

power_means <- list(
  colMeans(power1[,1:len], na.rm = TRUE),
  colMeans(power2[,1:len], na.rm = TRUE),
  colMeans(power3[,1:len], na.rm = TRUE),
  colMeans(power4[,1:len], na.rm = TRUE),
  colMeans(power5[,1:len], na.rm = TRUE),
  colMeans(power6[,1:len], na.rm = TRUE),
  colMeans(power7[,1:len], na.rm = TRUE)
)

fdr_means <- list(
  colMeans(fdr1[,1:len], na.rm = TRUE),
  colMeans(fdr2[,1:len], na.rm = TRUE),
  colMeans(fdr3[,1:len], na.rm = TRUE),
  colMeans(fdr4[,1:len], na.rm = TRUE),
  colMeans(fdr5[,1:len], na.rm = TRUE),
  colMeans(fdr6[,1:len], na.rm = TRUE),
  colMeans(fdr7[,1:len], na.rm = TRUE)
)

df_power <- do.call(rbind, lapply(1:7, function(i) {
  data.frame(Method = methods[i],
             q_value = fdr_seq,
             Value = power_means[[i]],
             Metric = "Power")
}))

df_fdr <- do.call(rbind, lapply(1:7, function(i) {
  data.frame(Method = methods[i],
             q_value = fdr_seq,
             Value = fdr_means[[i]],
             Metric = "FDR")
}))

df <- rbind(df_power, df_fdr)


df$Metric <- factor(df$Metric, levels = c("Power", "FDR"))

df <- df[df$Method %in% c("Knockoff", "AnnoKn-simple", "AnnoKn", 
                          "GhostKnockoff", "AnnoGK", "AdaKn (RF)"),]


ggplot(df, aes(x = q_value, y = Value, 
               color = Method, shape = Method, linetype = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~Metric, ncol = 2, scales = "free_y") +
  geom_line(
    data = data.frame(
      q_value = fdr_seq, 
      Value = fdr_seq, 
      Metric = factor("FDR", levels = c("Power","FDR"))
    ),
    aes(x = q_value, y = Value),
    inherit.aes = FALSE, 
    color = "black", linetype = "dashed"
  ) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_minimal(base_size = 14) +
  labs(
    x = "q-value threshold", 
    y = "Mean Value", 
    title = paste0("n_", N.effect, "_h2_", h2e, "_alphae_", alphae, "_p1_", p1, "_regu_", regu)
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )






####################### for multiple h^2 #################






N.effect = 10000 # sample size
n1 = 1000 # sample size used for clustering
# p = 200 # number of variable
h2e = 0.02 # heritability 0.02, 0.05, 0.1, 0.2
alphae = 0.3 # proportion of causal SNPs. 0.05, 0.1, 0.2
# rho = 0.5 # maximum correlation between SNPs
M = 1 # number of knockoff copies
chrid = '1'
regu = 0.2
width = "2MB"
if(width == '2MB'){
  p1 = 100
}else{
  p1 = 200
}

df_all <- c()

for(h2e in c(0.02, 0.05, 0.1)){
  
  power1 <- power2 <- power3 <- power4 <- power5 <- power6 <- power7 <- numeric()
  fdr1 <- fdr2 <- fdr3 <- fdr4 <- fdr5 <- fdr6 <- fdr7<- numeric()
  lambda <- numeric()
  
  for (iter in 1:100) {
    
    path = paste0(repro_path("true_genotype_paper_results", "result_"), width, "_insample_iter_", iter, "_n_", N.effect, "_h2e_", h2e, "_alphae_", alphae, "_p1_", p1, "_M_", M, "_chr_", chrid, ".RData")
    load(path)
    
    power1 <- rbind(power1, result$power1)
    power2 <- rbind(power2, result$power2)
    power3 <- rbind(power3, result$power3)
    power4 <- rbind(power4, result$power4)
    power5 <- rbind(power5, result$power5)
    power6 <- rbind(power6, result$power6)
    power7 <- rbind(power7, result$power7)
    
    fdr1 <- rbind(fdr1, result$fdr1)
    fdr2 <- rbind(fdr2, result$fdr2)
    fdr3 <- rbind(fdr3, result$fdr3)
    fdr4 <- rbind(fdr4, result$fdr4)
    fdr5 <- rbind(fdr5, result$fdr5)
    fdr6 <- rbind(fdr6, result$fdr6)
    fdr7 <- rbind(fdr7, result$fdr7)
    
    lambda = append(lambda, result$lambda)
  }
  
  
  library(ggplot2)
  
  methods <- c("Knockoff", "AnnoKn-lite", "AnnoKn", 
               "GhostKnockoff", "AnnoGK-lite", "AnnoGK", "AdaKn (RF)")
  
  fdr_seq <- seq(0.4, 0.1, by = -0.05)
  len = length(fdr_seq)
  
  power_means <- list(
    colMeans(power1[,1:len], na.rm = TRUE),
    colMeans(power2[,1:len], na.rm = TRUE),
    colMeans(power3[,1:len], na.rm = TRUE),
    colMeans(power4[,1:len], na.rm = TRUE),
    colMeans(power5[,1:len], na.rm = TRUE),
    colMeans(power6[,1:len], na.rm = TRUE),
    colMeans(power7[,1:len], na.rm = TRUE)
  )
  
  fdr_means <- list(
    colMeans(fdr1[,1:len], na.rm = TRUE),
    colMeans(fdr2[,1:len], na.rm = TRUE),
    colMeans(fdr3[,1:len], na.rm = TRUE),
    colMeans(fdr4[,1:len], na.rm = TRUE),
    colMeans(fdr5[,1:len], na.rm = TRUE),
    colMeans(fdr6[,1:len], na.rm = TRUE),
    colMeans(fdr7[,1:len], na.rm = TRUE)
  )
  
  df_power <- do.call(rbind, lapply(1:7, function(i) {
    data.frame(Method = methods[i],
               q_value = fdr_seq,
               Value = power_means[[i]],
               Metric = "Power")
  }))
  
  df_fdr <- do.call(rbind, lapply(1:7, function(i) {
    data.frame(Method = methods[i],
               q_value = fdr_seq,
               Value = fdr_means[[i]],
               Metric = "FDR")
  }))
  
  df <- rbind(df_power, df_fdr)
  
  
  df$Metric <- factor(df$Metric, levels = c("Power", "FDR"))
  
  df <- df[df$Method %in% c("Knockoff", "AnnoKn-lite", "AnnoKn", 
                            "GhostKnockoff", "AnnoGK"),]
  df$heritability = paste0("h\u00B2 = ", h2e)
  
  df_all <- rbind(df_all, df)
  
}




# 
# colors <- c("Knockoff" = "#E69F00",
#             "AnnoKn" = "#FF0000",
#             "GhostKnockoff" = "#7570b3",
#             "AnnoGK" = "#e7298a",
#             "AnnoKn-lite" = "#66a61e",
#             "AnnoGK-lite" = "#e6ab02",
#             "AdaKn (RF)" = "#1b9e77")
# 
# shapes <- c("Knockoff" = 16,
#             "AnnoKn" = 8,
#             "GhostKnockoff" = 15,
#             "AnnoGK" = 3,
#             "AnnoKn-lite" = 5,
#             "AnnoGK-lite" = 4,
#             "AdaKn (RF)" = 17)
# 
# linetypes <- c(
#   "Knockoff" = "solid",
#   "AnnoKn" = "dashed",
#   "GhostKnockoff" = "dotdash",
#   "AnnoGK" = "twodash",
#   "AnnoKn-lite" = "dotted",
#   "AnnoGK-lite" = "longdash",
#   "AdaKn (RF)" = "solid"
# )


colors <- c("Knockoff" = "#e6ab02",
            "AnnoKn" = "#F8766D",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#4DAF4A",
            'AdaKn (RF)' = '#e7298a',
            "AnnoKn-lite" = "#619CFF",
            "AnnoGK-lite" = "#D55E00")

shapes <- c("Knockoff" = 16,
            "AnnoKn" = 8,
            "GhostKnockoff" = 15,
            "AnnoGK" = 3,
            'AdaKn (RF)' = 17,
            "AnnoKn-lite" = 5,
            "AnnoGK-lite" = 4)


linetypes <- c("Knockoff" = "solid",
               "AnnoKn" = "dashed",
               "GhostKnockoff" = "23",
               "AdaKn (RF)" = "22",
               "AnnoGK" = "21",
               "AnnoKn-lite" = "dotted",
               "AnnoGK-lite" = "24")




plot_data_with_limits <- df_all %>%
  mutate(ymax = ifelse(Metric == "FDR", 0.5, 0.7))

g <- ggplot(df_all, aes(x = q_value, y = Value, color = Method, 
                        shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  facet_grid(Metric ~ heritability, switch = "y", scales = "free_y") +  # 允许不同行不同y轴
  theme_minimal() +
  labs(x = "Target FDR", y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
  )

fdr_data <- subset(df_all, Metric == "FDR")
g <- g + geom_line(
  data = fdr_data,
  aes(x = q_value, y = q_value),
  inherit.aes = FALSE,
  linetype = "solid",
  color = "black",
  linewidth = 1
)

#g <- g + geom_blank(
#  data = plot_data_with_limits,
#  aes(y = ymax)
#)

print(g)


ggsave(repro_path("plots", "AnnoGK_simu_true_genotype_insample.pdf"),
       g, 
       width = 9, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)

