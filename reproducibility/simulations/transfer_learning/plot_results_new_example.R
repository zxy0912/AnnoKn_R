if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)

Sys.setenv(TZ = 'America/New_York')
library(tidyverse)
library(gridExtra)

N = 1
type = 'single'
error = 'corrected'
N_method = 11
FDPs = NULL
POWERs = NULL
Methods = NULL
lambda_l = c()
for (arg in (1:300)){
  filename = paste(repro_path("transfer_results", "new_output_"),type, error, arg,".Rdata", sep = "")
  if(file.exists(filename)){
    load(filename)
    FDPs = c(FDPs, fdps)
    POWERs = c(POWERs, powers)
    Methods = c(Methods, methods)
    lambda_l = rbind(lambda_l, unlist(lambda_s))
  }
}
results = data.frame(power = POWERs, fdp = FDPs, method = Methods)

method.values <- c("pool", "vanilla", "linear_com_0.1", "linear_com_0.2", "linear_com_0.3", "linear_com_0.4", "linear_com_0.5", "gam", "weighted_lasso", "AnnoKn")
method.labels <- c("Heuristic (pool)", "Vanilla", "Linear (0.1)", "Linear (0.2)", "Linear (0.3)", "Linear (0.4)", "Linear (0.5)", "Transfer (gam)", "Weighted-lasso", "AnnoKn")

# 2.
plot_data <- results %>%

  filter(method %in% method.values) %>%
  mutate(Method = factor(method, levels = method.values, labels = method.labels)) %>%
  rename(Power = power, FDR = fdp) %>%
  pivot_longer(cols = c(Power, FDR), names_to = "Metric", values_to = "Value")

# 3.
df.nominal <- data.frame(Metric = "FDR", Value = 0.1)

# 4.
p_box <- ggplot(plot_data, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  geom_hline(data = df.nominal, aes(yintercept = Value), linetype = "dashed", color = "red", size = 0.8) +
  scale_fill_manual(values = c(
    "Heuristic (pool)" = "#377EB8",
    "Vanilla" = "#4DAF4A",
    "Linear (0.1)" = "#984EA3",
    "Linear (0.2)" = "#984EA3",
    "Linear (0.3)" = "#984EA3",
    "Linear (0.4)" = "#984EA3",
    "Linear (0.5)" = "#984EA3",
    "Transfer (gam)" = "#FF7F00",
    "Weighted-lasso" = "#A65628",
    "AnnoKn" = "#E41A1C"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

print(p_box)








########## combining two AnnoKn methods:


Sys.setenv(TZ = 'America/New_York')
library(tidyverse)
library(gridExtra)

N = 1
error = 'corrected'
N_method = 11

# ==========================================
# 1.
# ==========================================
FDPs_all = NULL
POWERs_all = NULL
Methods_all = NULL
lambda_l = c()

for (arg in (1:300)){
  filename = paste(repro_path("transfer_results", "new_output_"),'all', error, arg,".Rdata", sep = "")
  if(file.exists(filename)){
    load(filename)
    FDPs_all = c(FDPs_all, fdps)
    POWERs_all = c(POWERs_all, powers)
    methods[methods == "AnnoKn"] <- "AnnoKn_all"
    Methods_all = c(Methods_all, methods)
    lambda_l = rbind(lambda_l, unlist(lambda_s))
  }
}
df_all = data.frame(power = POWERs_all, fdp = FDPs_all, method = Methods_all)

# ==========================================
# 2.
# ==========================================
FDPs_single = NULL
POWERs_single = NULL
Methods_single = NULL

annokn_idx_all = c()

for (arg in (1:300)){
  filename = paste(repro_path("transfer_results", "new_output_"),'single', error, arg,".Rdata", sep = "")
  if(file.exists(filename)){
    load(filename)

    annokn_idx = which(methods == "AnnoKn")
    annokn_idx_all = append(annokn_idx_all, annokn_idx)
    FDPs_single = c(FDPs_single, fdps[annokn_idx])
    POWERs_single = c(POWERs_single, powers[annokn_idx])
    Methods_single = c(Methods_single, rep("AnnoKn_single", length(annokn_idx)))
  }
}
df_single = data.frame(power = POWERs_single, fdp = FDPs_single, method = Methods_single)

# ==========================================
# 3.
# ==========================================
results = rbind(df_all, df_single)

method.values <- c("pool", "vanilla", "linear_com_0.1", "linear_com_0.2", "linear_com_0.3", "linear_com_0.4", "linear_com_0.5", "gam", "weighted_lasso", "AnnoKn_single", "AnnoKn_all")
method.labels <- c("Heuristic (pool)", "Vanilla", "Linear (0.1)", "Linear (0.2)", "Linear (0.3)", "Linear (0.4)", "Linear (0.5)", "Transfer (gam)", "Weighted-lasso", "AnnoKn-1D", "AnnoKn-2D")

plot_data <- results %>%
  filter(method %in% method.values) %>%
  mutate(Method = factor(method, levels = method.values, labels = method.labels)) %>%
  rename(Power = power, FDR = fdp) %>%
  pivot_longer(cols = c(Power, FDR), names_to = "Metric", values_to = "Value")

df.nominal <- data.frame(Metric = "FDR", Value = 0.1)

# ==========================================
# 4.
# ==========================================
p_box <- ggplot(plot_data, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  geom_hline(data = df.nominal, aes(yintercept = Value), linetype = "dashed", color = "red", size = 0.8) +

  scale_fill_manual(values = c(
    "Heuristic (pool)" = "#377EB8",
    "Vanilla" = "#4DAF4A",
    "Linear (0.1)" = "#984EA3",
    "Linear (0.2)" = "#984EA3",
    "Linear (0.3)" = "#984EA3",
    "Linear (0.4)" = "#984EA3",
    "Linear (0.5)" = "#984EA3",
    "Transfer (gam)" = "#FF7F00",
    "Weighted-lasso" = "#A65628",
    "AnnoKn-1D" = "#FB9A99",
    "AnnoKn-2D" = "#E41A1C"
  )) +

  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

print(p_box)




ggsave(
  filename = repro_path("transfer_plots", "diff_overlap.pdf"),
  plot = p_box,
  width = 10,
  height = 6,
  device = "pdf",
  create.dir = TRUE
)



############

library(tidyverse)

# 1.
df_lambda <- as.data.frame(lambda_l)
colnames(df_lambda) <- c("External Env 1 (20% Overlap)", "External Env 2 (80% Overlap)")

# 2.
plot_data_lambda <- df_lambda %>%
  pivot_longer(cols = everything(), names_to = "Environment", values_to = "Lambda")

# 3.
p_lambda_hist <- ggplot(plot_data_lambda, aes(x = Lambda, fill = Environment)) +
  geom_histogram(position = "identity", alpha = 0.6, color = "black", bins = 50) +
  scale_fill_manual(values = c(
    "External Env 1 (20% Overlap)" = "#377EB8",
    "External Env 2 (80% Overlap)" = "#E41A1C"
  )) +
  theme_bw() +
  labs(
    x = expression(lambda[l]),
    y = "Count"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

print(p_lambda_hist)


ggsave(
  filename = repro_path("transfer_plots", "lambda.pdf"),
  plot = p_lambda_hist,
  width = 8,
  height = 5,
  device = "pdf",
  create.dir = TRUE
)
