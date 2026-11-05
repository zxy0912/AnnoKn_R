if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)
Sys.setenv(TZ = 'America/New_York')
library(tidyverse)
library(gridExtra)

N = 10
N_method = 11
type = 'single'
error = 'corrected'
FDPs = NULL
POWERs = NULL
overlaps = NULL
Methods = NULL
for (arg in (1:300)){
  filename = paste(repro_path("transfer_results", "output_"),type, error, arg,".Rdata", sep = "")
  # filename = paste(repro_path("transfer_results", "output"), arg,".Rdata", sep = "")

  if(file.exists(filename)){
    load(filename)
    overlap = 20*((arg-1)%%6)
    FDPs = c(FDPs, fdps)
    POWERs = c(POWERs, powers)
    overlaps = c(overlaps, rep(overlap,N*N_method))
    Methods = c(Methods, methods)
  }
}
results = data.frame(overlap = overlaps, power = POWERs, fdp = FDPs, method = Methods)

df.nominal <- tibble(key="FDR", value=0.1)

method.values <- c("pool", "vanilla", "adaptive", paste("linear_com", 0.1*(1:5), sep = "_"), "oracle", "gam", "weighted_lasso", "AnnoKn")
method.labels <- c("Heuristic (pool)", "Vanilla", "Adaptive (EM)",
                   paste("Transfer - linear_com", 0.1*(1:5), sep = "_"), #rep("Transfer - linearly re-ordered",5),
                   "Transfer - linearly re-ordered (oracle)" , "Transfer - adaptive (gam)", "Transfer - weighted-lasso", "AnnoKn")
#method.labels <- c("Heuristic (pool)", "Vanilla knockoffs", "Adaptive knockoffs", paste("Linear combination statistics, theta = ", 0.1*(1:5), sep = ""),"Weighted Lasso Statistics")

## find the oracle method
# NN = 300
# power_matrix = cbind(POWERs[10*(0:(NN-1)) + 5],POWERs[10*(0:(NN-1)) + 6],POWERs[10*(0:(NN-1)) + 7], POWERs[10*(0:(NN-1)) + 8], POWERs[10*(0:(NN-1)) + 9])
# which_theta = apply(power_matrix, 1, which.max)
# oracle_index = 10*(0:(NN-1)) + which_theta + 4
# results_oracle = data.frame(overlap = overlaps[oracle_index], power = POWERs[oracle_index], fdp = FDPs[oracle_index], method = rep("oracle", NN))
# results = rbind(results, results_oracle)


linear_summary = results %>%
  filter(method == "linear_com_0.1"|
           method == "linear_com_0.2"|
           method == "linear_com_0.3"|
           method == "linear_com_0.4"|
           method == "linear_com_0.5") %>%
  mutate(Power = power, FDR = fdp) %>%
  group_by(overlap, method) %>%
  summarise(FDR = mean(FDR), Power = mean(Power)) %>%
  group_by(overlap) %>%
  summarise(idx = which.max(Power), Power=Power[idx], FDR = FDR[idx]) %>%
  select(-idx) %>%
  gather(Power, FDR, key="key", value="value") %>%
  mutate(method = "oracle")

others_summary = results %>%
  filter(method != "linear_com_0.1") %>%
  filter(method != "linear_com_0.2") %>%
  filter(method != "linear_com_0.3") %>%
  filter(method != "linear_com_0.4") %>%
  filter(method != "linear_com_0.5") %>%
  filter(method != "adaptive") %>%
  mutate(Power = power, FDR = fdp) %>%
  gather(Power, FDR, key="key", value="value") %>%
  group_by(overlap, method, key) %>%
  summarise(num=n(), value.se=2*sd(value)/sqrt(num), value=mean(value, na.rm = TRUE))

full_summary = full_join(linear_summary, others_summary)


p1 = full_summary %>%
  mutate(Method = factor(method, method.values, method.labels)) %>%
  #mutate(value = ifelse(is.na(value), 0, value)) %>%
  ggplot(aes(x=overlap, y=value, color=Method, shape=Method)) +
  #ggplot(aes(x=overlap, y=value, color=Method)) +
  geom_point(alpha=0.8) +
  geom_line(alpha=0.8) +
  #geom_errorbar(aes(ymin=value-value.se, ymax=value+value.se), width=0.2) +
  geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
  scale_color_manual(values = c(
    "Heuristic (pool)" = "#377EB8",
    "Vanilla" = "#4DAF4A",
    "Transfer - linearly re-ordered (oracle)" = "#984EA3",
    "Transfer - adaptive (gam)" = "#FF7F00",
    "Transfer - weighted-lasso" = "#A65628",
    "AnnoKn" = "#E41A1C"
  )) +
  #scale_shape_manual(values=shape.scale) +
  facet_wrap(.~key, nrow=1, scales="free_y") +
  xlab("Overlap (in %)") +
  ylab(NULL) +
  ylim(0,1) +
  theme_bw() +
  theme(#strip.background = element_blank(),
    #strip.placement = "outside",
    legend.position = "right",
    plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10)) +
  guides(col = guide_legend(nrow = 7)) #+
#ggtitle("(a) Setting 1 (high invariance)")

p1
# p1 %>% ggsave(file="experiment_transfer.pdf", width=7.8, height=2.7, units="in")




###########

FDPs = NULL
POWERs = NULL
overlaps = NULL
Methods = NULL
for (arg in (1:300)){
  filename = paste(repro_path("transfer_results", "output_"),type, error, arg,".Rdata", sep = "")
  if(file.exists(filename)){
    load(filename)
    overlap = 20*((arg-1)%%6)
    FDPs = c(FDPs, fdps)
    POWERs = c(POWERs, powers)
    overlaps = c(overlaps, rep(overlap,N*N_method))
    Methods = c(Methods, methods)
  }
}
results = data.frame(overlap = overlaps, power = POWERs, fdp = FDPs, method = Methods)

df.nominal <- tibble(key="FDR", value=0.1)

method.values <- c("pool", "vanilla", "linear_com_0.2", "gam", "weighted_lasso", "AnnoKn")
method.labels <- c("Heuristic (pool)", "Vanilla", "Transfer - linear_theta_0.2",
                   "Transfer - adaptive (gam)", "Transfer - weighted-lasso", "AnnoKn")

## find the oracle method
# NN = 300
# power_matrix = cbind(POWERs[10*(0:(NN-1)) + 5],POWERs[10*(0:(NN-1)) + 6],POWERs[10*(0:(NN-1)) + 7], POWERs[10*(0:(NN-1)) + 8], POWERs[10*(0:(NN-1)) + 9])
# which_theta = apply(power_matrix, 1, which.max)
# oracle_index = 10*(0:(NN-1)) + which_theta + 4
# results_oracle = data.frame(overlap = overlaps[oracle_index], power = POWERs[oracle_index], fdp = FDPs[oracle_index], method = rep("oracle", NN))
# results = rbind(results, results_oracle)



full_summary = results %>%
  filter(method != "linear_com_0.1") %>%
  filter(method != "linear_com_0.3") %>%
  filter(method != "linear_com_0.4") %>%
  filter(method != "linear_com_0.5") %>%
  filter(method != "adaptive") %>%
  mutate(Power = power, FDR = fdp) %>%
  gather(Power, FDR, key="key", value="value") %>%
  group_by(overlap, method, key) %>%
  summarise(num=n(), value.se=2*sd(value)/sqrt(num), value=mean(value, na.rm = TRUE))



p2 = full_summary %>%
  mutate(Method = factor(method, method.values, method.labels)) %>%
  #mutate(value = ifelse(is.na(value), 0, value)) %>%
  ggplot(aes(x=overlap, y=value, color=Method, shape=Method)) +
  #ggplot(aes(x=overlap, y=value, color=Method)) +
  geom_point(alpha=0.8) +
  geom_line(alpha=0.8) +
  #geom_errorbar(aes(ymin=value-value.se, ymax=value+value.se), width=0.2) +
  geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
  scale_color_manual(values = c(
    "Heuristic (pool)" = "#377EB8",                            # 蓝色
    "Vanilla" = "#4DAF4A",                                     # 绿色
    "Transfer - linear_theta_0.2" = "#984EA3",     # 紫色
    "Transfer - adaptive (gam)" = "#FF7F00",                   # 橙色
    "Transfer - weighted-lasso" = "#A65628",                   # 棕色 (替代了原本的黄色)
    "AnnoKn" = "#E41A1C"                                       # 红色 (作为你们的方法，用最醒目的红色突出)
  )) +
  #scale_shape_manual(values=shape.scale) +
  facet_wrap(.~key, nrow=1, scales="free_y") +
  xlab("Overlap (in %)") +
  ylab(NULL) +
  ylim(0,1) +
  theme_bw() +
  theme(#strip.background = element_blank(),
    #strip.placement = "outside",
    legend.position = "right",
    plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10)) +
  guides(col = guide_legend(nrow = 7)) #+
#ggtitle("(a) Setting 1 (high invariance)")

p2
# p1 %>% ggsave(file="experiment_transfer.pdf", width=7.8, height=2.7, units="in")

ggsave(
  filename = repro_path("transfer_plots", "benchmark.pdf"),
  plot = p2,
  width = 8,
  height = 3.5,
  device = "pdf",
  create.dir = TRUE
)
