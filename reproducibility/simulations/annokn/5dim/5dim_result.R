# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)





library(ggplot2)

iteration = 120
k = 10

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)
num_noise = 4

############################ 1 dimension annotation ############################## 

power1 = matrix(0, nrow = iteration, ncol = len)
power2 = matrix(0, nrow = iteration, ncol = len)
power3 = matrix(0, nrow = iteration, ncol = len)
power4 = matrix(0, nrow = iteration, ncol = len)
power5 = matrix(0, nrow = iteration, ncol = len)
power6 = matrix(0, nrow = iteration, ncol = len)
power7 = matrix(0, nrow = iteration, ncol = len)
power8 = matrix(0, nrow = iteration, ncol = len)
power9 = matrix(0, nrow = iteration, ncol = len)

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)
fdr9 = matrix(0, nrow = iteration, ncol = len)

lambdas_annokn <- matrix(0, nrow = iteration, ncol = num_noise+1)
lambdas_annokn_lite <- matrix(0, nrow = iteration, ncol = num_noise+1)

causal_index <- matrix(0, nrow = iteration, ncol = 1)


s = 0
removelist = c()

n = 1000
p = 900
k = 150
d = 2
amp = 3.5
binary = 'binary'

for (i in 1:iteration) {
  path = paste0(repro_path("annokn_5d_results", "result_randomposition_1_p_"), num_noise+1, "_",binary, "_amp_",amp, "_", i,".RData")
  if (!file.exists(path)){
    print(paste("iteration",i,"doesn't exist"))
    removelist = append(removelist, i)
    next
  }
  load(path) 
  
  power1[i,] <- result$power1
  power2[i,] <- result$power2
  power3[i,] <- result$power3
  power4[i,] <- result$power4
  power5[i,] <- result$power5
  power6[i,] <- result$power6
  power7[i,] <- result$power7
  power8[i,] <- result$power8
  power9[i,] <- result$power9
  
  
  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8
  fdr9[i,]   <- result$fdr9
  
  lambdas_annokn[i,] <- result$lambda_vals_annakn
  lambdas_annokn_lite[i,] <- result$lambda_vals_annakn_lite
  causal_index[i,] <- result$causal_anno_index
  
}

if(length(removelist) > 0){
  power1 = power1[-removelist,]
  power2 = power2[-removelist,]
  power3 = power3[-removelist,]
  power4 = power4[-removelist,]
  power5 = power5[-removelist,]
  power6 = power6[-removelist,]
  power7 = power7[-removelist,]
  power8 = power8[-removelist,]
  power9 = power9[-removelist,]
  fdr1 = fdr1[-removelist,]
  fdr2 = fdr2[-removelist,]
  fdr3 = fdr3[-removelist,]
  fdr4 = fdr4[-removelist,]
  fdr5 = fdr5[-removelist,]
  fdr6 = fdr6[-removelist,]
  fdr7 = fdr7[-removelist,]
  fdr8 = fdr8[-removelist,]
  fdr9 = fdr9[-removelist,]
  lambdas_annokn <- lambdas_annokn[-removelist,]
  lambdas_annokn_lite <- lambdas_annokn_lite[-removelist,]
  causal_index <- as.matrix(causal_index[-removelist,], ncol = 1)
}

number_totel <- 100 - length(removelist)

power <- c(
  colMeans(power1[1:100,], na.rm = TRUE),
  colMeans(power2[1:100,], na.rm = TRUE),
  colMeans(power3[1:100,], na.rm = TRUE),
  colMeans(power8[1:100,], na.rm = TRUE),  
  colMeans(power9[1:100,], na.rm = TRUE),
  colMeans(power4[1:100,], na.rm = TRUE),  
  colMeans(power5[1:100,], na.rm = TRUE),
  colMeans(power6[1:100,], na.rm = TRUE),  
  colMeans(power7[1:100,], na.rm = TRUE)
)


fdr <- c(
  colMeans(fdr1[1:100,], na.rm = TRUE),
  colMeans(fdr2[1:100,], na.rm = TRUE),
  colMeans(fdr3[1:100,], na.rm = TRUE),
  colMeans(fdr8[1:100,], na.rm = TRUE),
  colMeans(fdr9[1:100,], na.rm = TRUE),
  colMeans(fdr4[1:100,], na.rm = TRUE),
  colMeans(fdr5[1:100,], na.rm = TRUE),
  colMeans(fdr6[1:100,], na.rm = TRUE),
  colMeans(fdr7[1:100,], na.rm = TRUE)
)

lambdas_annokn <- lambdas_annokn[1:100,]
lambdas_annokn_lite <- lambdas_annokn_lite[1:100,]
causal_index <- causal_index[1:100,]


name_method <- c('Knockoff', 
                 'AdaKn (EM) with informative annotation',
                 'AdaKn (EM) with all annotation', 
                 'AdaKn (RF) with informative annotation',
                 'AdaKn (RF) with all annotation',
                 'AnnoKn with informative annotation', 
                 'AnnoKn with all annotation', 
                 'AnnoKn-lite with informative annotation',
                 'AnnoKn-lite with all annotation')

methods = name_method

method <- rep(name_method, each = len)
plot_data <- data.frame(alpha = rep(alphalist, length(name_method)),
                        method = method,
                        power = power,
                        fdr = fdr)


plot_data <- plot_data[plot_data$method %in% c('Knockoff',
                                               'AnnoKn with informative annotation',
                                               'AnnoKn with all annotation',
                                               'AnnoKn-lite with informative annotation',
                                               'AnnoKn-lite with all annotation',
                                               'AdaKn (EM) with informative annotation',
                                               'AdaKn (EM) with all annotation'),]

colors <- c('Knockoff' = '#e6ab02', 
            'AdaKn (EM) with informative annotation' = '#C77CFF',
            'AdaKn (EM) with all annotation' = '#C77CFF',
            'AnnoKn with informative annotation'= '#F8766D', 
            'AnnoKn with all annotation'= '#F8766D', 
            'AnnoKn-lite with informative annotation' = '#619CFF',
            'AnnoKn-lite with all annotation' = '#619CFF', 
            'AdaKn (RF) with informative annotation' = '#e7298a',
            'AdaKn (RF) with all annotation' = '#e7298a')

shapes <- c("Knockoff" = 16,
            'AdaKn (EM) with informative annotation' = 7,
            'AdaKn (EM) with all annotation' = 7,
            'AnnoKn with informative annotation'= 8, 
            'AnnoKn with all annotation'= 8, 
            'AnnoKn-lite with informative annotation' = 5,
            'AnnoKn-lite with all annotation' = 5, 
            'AdaKn (RF) with informative annotation' = 17,
            'AdaKn (RF) with all annotation' = 17)

linetypes <- c("Knockoff" = "solid", # 1
               'AdaKn (EM) with informative annotation' = "solid",
               'AdaKn (EM) with all annotation' = "dashed",
               'AnnoKn with informative annotation'= "solid", 
               'AnnoKn with all annotation'= "dashed", 
               'AnnoKn-lite with informative annotation' = "solid",
               'AnnoKn-lite with all annotation' = "dashed", 
               'AdaKn (RF) with informative annotation' = "solid",
               'AdaKn (RF) with all annotation' = "dashed") 


get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1_with_legend <- ggplot(plot_data, aes(x = alpha, y = power, color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, breaks = methods) +
  scale_shape_manual(values = shapes, breaks = methods) +
  scale_linetype_manual(values = linetypes, breaks = methods) +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 10)) 

my_legend <- get_legend(p1_with_legend)

p1 <- ggplot(plot_data, aes(x = alpha, y = power, color = method, shape = method, linetype = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Target FDR", y = "Power") +
  scale_color_manual(values = colors, breaks = methods) +
  scale_shape_manual(values = shapes, breaks = methods) +
  scale_linetype_manual(values = linetypes, breaks = methods) +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- ggplot(plot_data, aes(x = alpha, y = fdr, color = method, shape = method, linetype = method)) +
  geom_line() +
  geom_point() +
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_segment(
    aes(x = 0.05, y = 0.05, xend = 0.3, yend = 0.3), 
    linetype = "solid", 
    color = "black", 
    linewidth = 1
  ) +
  labs(x = "Target FDR", y = "FDR") +
  scale_color_manual(values = colors, breaks = methods) +
  scale_shape_manual(values = shapes, breaks = methods) +
  scale_linetype_manual(values = linetypes, breaks = methods) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0.05, 0.3), ylim = c(0, 0.3))

library(gridExtra)
g <- arrangeGrob(p1, p2, my_legend, ncol = 3, widths = c(1, 1, 0.8))


library(grid)
grid.newpage()
grid.draw(g)


ggsave(repro_path("plots", "AnnoKn_5dimen_binary.pdf"),
       g, 
       width = 11, height = 4, units = "in", 
       bg = "white", device = cairo_pdf)





lambda_inf <- mapply(function(row, col) lambdas_annokn[row, col],
                     row = 1:nrow(lambdas_annokn),
                     col = causal_index)

# Step 2: 提取 non-informative 参数（该行除 causal_index 那一列外）
lambda_noninf <- mapply(function(row, col) {
  lambdas_annokn[row, -col]
}, row = 1:nrow(lambdas_annokn),
col = causal_index, SIMPLIFY = FALSE)

# 转成向量
lambda_inf_vec    <- as.numeric(lambda_inf)
lambda_noninf_vec <- as.numeric(unlist(lambda_noninf))

# Step 3: 组合数据框
df <- rbind(
  data.frame(value = lambda_noninf_vec, group = "non-informative"),
  data.frame(value = lambda_inf_vec,    group = "informative")
)

# Step 4: 绘制分布图（直方图，自适应范围）
g <- ggplot(na.omit(df), aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5,
                 bins = 60) +  # 自动等宽分50个bin，可根据数据调整
  scale_fill_manual(values = c("non-informative" = "#E69F00",  # 橙
                               "informative"     = "#0072B2")) +  # 蓝
  labs(x = expression(lambda[l]), y = "Frequency", fill = "Group") +
  theme_minimal(base_size = 13)

ggsave(repro_path("plots", "Distribution_5dimen_continuous.pdf"),
       g, 
       width = 10, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)


