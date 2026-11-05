if (!exists("repro_path", mode = "function")) stop("Use Rscript reproducibility/run.R with this script path.", call. = FALSE)


library(ggplot2)

iteration = 150

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)

############################ 1 dimension annotation ##############################

power1 = matrix(0, nrow = iteration, ncol = len)
power2 = matrix(0, nrow = iteration, ncol = len)
power3 = matrix(0, nrow = iteration, ncol = len)
power4 = matrix(0, nrow = iteration, ncol = len)
power5 = matrix(0, nrow = iteration, ncol = len)
power6 = matrix(0, nrow = iteration, ncol = len)
power7 = matrix(0, nrow = iteration, ncol = len)
power8 = matrix(0, nrow = iteration, ncol = len)

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)

time1 <- time2 <- time3 <- time4 <- time5 <- time6 <- time7 <- time8 <- numeric()


s = 0
removelist = c()

n = 1000
p = 900
k = 150
d = 2
amp = 3.5
binary = ''
theta = 1

for (i in 1:iteration) {
  path = paste0(repro_path("annotation_dependence_results", "result_indep_"), binary, "_theta_",theta, "_amp_",amp, "_", i,".RData")
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


  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8

  time1 <- append(time1, result$time1)
  time2 <- append(time2, result$time2)
  time3 <- append(time3, result$time3)
  time4 <- append(time4, result$time4)
  time5 <- append(time5, result$time5)
  time6 <- append(time6, result$time6)
  time7 <- append(time7, result$time7)
  time8 <- append(time8, result$time8)

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
  fdr1 = fdr1[-removelist,]
  fdr2 = fdr2[-removelist,]
  fdr3 = fdr3[-removelist,]
  fdr4 = fdr4[-removelist,]
  fdr5 = fdr5[-removelist,]
  fdr6 = fdr6[-removelist,]
  fdr7 = fdr7[-removelist,]
  fdr8 = fdr8[-removelist,]
}


power <- c(
  colMeans(power1[1:100,], na.rm = TRUE),
  colMeans(power2[1:100,], na.rm = TRUE),
  colMeans(power3[1:100,], na.rm = TRUE),
  colMeans(power4[1:100,], na.rm = TRUE),
  colMeans(power5[1:100,], na.rm = TRUE),
  colMeans(power6[1:100,], na.rm = TRUE),
  colMeans(power7[1:100,], na.rm = TRUE),
  colMeans(power8[1:100,], na.rm = TRUE)
)


fdr <- c(
  colMeans(fdr1[1:100,], na.rm = TRUE),
  colMeans(fdr2[1:100,], na.rm = TRUE),
  colMeans(fdr3[1:100,], na.rm = TRUE),
  colMeans(fdr4[1:100,], na.rm = TRUE),
  colMeans(fdr5[1:100,], na.rm = TRUE),
  colMeans(fdr6[1:100,], na.rm = TRUE),
  colMeans(fdr7[1:100,], na.rm = TRUE),
  colMeans(fdr8[1:100,], na.rm = TRUE)
)

time1 <- time1[1:100]
time2 <- time2[1:100]
time3 <- time3[1:100]
time4 <- time4[1:100]
time5 <- time5[1:100]
time6 <- time6[1:100]
time7 <- time7[1:100]
time8 <- time8[1:100]


name_method <- c('AdaPT',
                 'Knockoff',
                 'AdaKn (GLM)',
                 'AdaKn (GAM)',
                 'AdaKn (RF)',
                 'AdaKn (EM)',
                 'AnnoKn',
                 'AnnoKn-lite')

methods = name_method

method <- rep(name_method, each = len)
plot_data <- data.frame(alpha = rep(alphalist, 8),
                        method = method,
                        power = power,
                        fdr = fdr)

if(binary == ''){
  plot_data <- plot_data[plot_data$method %in% c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)',
                                                 'AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite'),]
  methods = c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)','AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite')
}

if(binary == 'binary'){
  plot_data <- plot_data[plot_data$method %in% c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)',
                                                 'AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite'),]
  methods = c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)','AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite')
}




if(binary == ''){
  plot_data <- plot_data[plot_data$method %in% c('Knockoff', 'AnnoKn','AnnoKn-lite'),]
  methods = c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)','AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite')
}

if(binary == 'binary'){
  plot_data <- plot_data[plot_data$method %in% c('Knockoff', 'AnnoKn','AnnoKn-lite'),]
  methods = c('Knockoff', 'AdaKn (GLM)', 'AdaKn (GAM)','AdaKn (RF)','AdaKn (EM)','AnnoKn','AnnoKn-lite')
}




colors <- c('AdaPT' = '#1b9e77',
            'Knockoff' = '#e6ab02',
            'AdaKn (GLM)' = '#d95f02',
            'AdaKn (GAM)' = '#00BFC4',
            'AdaKn (RF)' = '#e7298a',
            'AdaKn (EM)' = '#C77CFF',
            'AnnoKn' = '#F8766D',
            'AnnoKn-lite' = '#619CFF')

shapes <- c('AdaPT' = 11,
            "Knockoff" = 16,
            'AdaKn (GLM)' = 12,
            'AdaKn (GAM)' = 13,
            'AdaKn (RF)' = 17,
            'AdaKn (EM)' = 7,
            "AnnoKn" = 8,
            "AnnoKn-lite" = 5)

linetypes <- c('AdaPT' = "82",
               "Knockoff" = "solid", # 1
               'AdaKn (GLM)' = "41",
               'AdaKn (GAM)' = "14",
               "AdaKn (RF)" = "22",
               'AdaKn (EM)' = "81",
               "AnnoKn" = "dashed", # 2
               "AnnoKn-lite" = "dotted") # 3


get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1_with_legend <- ggplot(plot_data, aes(x = alpha, y = power, color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 1) +
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
  #geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", linewidth = 1) +
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
  coord_cartesian(xlim = c(0.05, 0.3), ylim = c(0, 0.4))

library(gridExtra)
g <- arrangeGrob(p1, p2, my_legend, ncol = 3, widths = c(1, 1, 0.6))


library(grid)
grid.newpage()
grid.draw(g)

ggsave(repro_path("annotation_dependence_plots", "AnnoKn_indep_setting_selected.pdf"),
# ggsave(repro_path("plots", "AnnoKn_AdaKn_setting_binary.pdf"),
       g,
       width = 10, height = 4, units = "in",
       bg = "white", device = cairo_pdf)



############# running time:
df <- data.frame(
  time = c(time1, time2, time3, time4, time5, time6, time7, time8),
  method = factor(rep(name_method, each = 100), levels = name_method)
)

g <- ggplot(df, aes(x = method, y = time, fill = method)) +
  geom_boxplot(outlier.size = 1.2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    name = "Time (seconds)",
    labels = function(x) paste0(x, " s")
  ) +
  labs(x = "Method") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

g

ggsave(repro_path("plots", "AnnoKn_AdaKn_time_continuous.pdf"),
       g,
       width = 6, height = 5, units = "in",
       bg = "white", device = cairo_pdf)
