# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)
library(ggplot2)

iteration = 100

alphalist <- seq(0.4, 0.1, by = -0.05)
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
power9 = matrix(0, nrow = iteration, ncol = len)
power10 = matrix(0, nrow = iteration, ncol = len)
power11 = matrix(0, nrow = iteration, ncol = len)

power_em = matrix(0, nrow = iteration, ncol = 1)
power_em_gk = matrix(0, nrow = iteration, ncol = 1)

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)
fdr9 = matrix(0, nrow = iteration, ncol = len)
fdr10 = matrix(0, nrow = iteration, ncol = len)
fdr11 = matrix(0, nrow = iteration, ncol = len)
fdr_em = matrix(0, nrow = iteration, ncol = 1)
fdr_em_gk =matrix(0, nrow = iteration, ncol = 1)


s = 0
removelist = c()

binary = ''
dimen = 'double'
heri = 0.1
N.effect = 5000 #
p = 300

for (i in 1:iteration) {
  path = paste0(repro_path("annogk_results", "heri_final_10_"),heri,"_n_",N.effect,"_p_",p, "_", i,".RData")
  if (!file.exists(path)){
    print(paste("iteration",i,"doesn't exist"))
    removelist = append(removelist, i)
    next
  }
  load(path) 
  s  = s+1
  
  power1[i,] <- result$power1
  power2[i,] <- result$power2
  power3[i,] <- result$power3
  power4[i,] <- result$power4
  power5[i,] <- result$power5
  power6[i,] <- result$power6
  power7[i,] <- result$power7
  power8[i,] <- result$power8
  power9[i,] <- result$power9
  power10[i,] <- result$power10
  power11[i,] <- result$power11
  
  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8
  fdr9[i,]   <- result$fdr9
  fdr10[i,]   <- result$fdr10
  fdr11[i,]   <- result$fdr11
  
}

# if(length(removelist) > 0){
#   power1 = power1[-removelist,]
#   power2 = power2[-removelist,]
#   power3 = power3[-removelist,]
#   power4 = power4[-removelist,]
#   power5 = power5[-removelist,]
#   power6 = power6[-removelist,]
#   power7 = power7[-removelist,]
#   power8 = power8[-removelist,]
#   fdr1 = fdr1[-removelist,]
#   fdr2 = fdr2[-removelist,]
#   fdr3 = fdr3[-removelist,]
#   fdr4 = fdr4[-removelist,]
#   fdr5 = fdr5[-removelist,]
#   fdr6 = fdr6[-removelist,]
#   fdr7 = fdr7[-removelist,]
#   fdr8 = fdr8[-removelist,]
# }


power <- c(
  colMeans(power1, na.rm = TRUE),
  colMeans(power2, na.rm = TRUE),
  colMeans(power3, na.rm = TRUE),
  colMeans(power4, na.rm = TRUE),  
  colMeans(power5, na.rm = TRUE),
  colMeans(power6, na.rm = TRUE),  
  colMeans(power7, na.rm = TRUE),
  colMeans(power8, na.rm = TRUE),
  colMeans(power9, na.rm = TRUE), 
  colMeans(power10, na.rm = TRUE),
  colMeans(power11, na.rm = TRUE) 
)


fdr <- c(
  colMeans(fdr1, na.rm = TRUE),
  colMeans(fdr2, na.rm = TRUE),
  colMeans(fdr3, na.rm = TRUE),
  colMeans(fdr4, na.rm = TRUE),
  colMeans(fdr5, na.rm = TRUE),
  colMeans(fdr6, na.rm = TRUE),
  colMeans(fdr7, na.rm = TRUE),
  colMeans(fdr8, na.rm = TRUE),
  colMeans(fdr9, na.rm = TRUE),
  colMeans(fdr10, na.rm = TRUE),
  colMeans(fdr11, na.rm = TRUE)
)


name_method <- c('Knockoff', 
                 'AnnoKn-simple', 
                 'AnnoKn', 
                 'GhostKnockoff', 
                 'AnnoGK-simple',
                 'AnnoGK',
                 'AnnoGK-dss', 
                 'GhostKnockoff M=5', 
                 'AnnoGK-simple M=5',
                 'AnnoGK M=5', 
                 'AnnoGK-dss M=5')

method <- rep(name_method, each = len)
plot_data_all <- data.frame(alpha = rep(alphalist, 11),
                        method = method,
                        power = power,
                        fdr = fdr)

plot_data <- plot_data_all[plot_data_all$method %in% c('Knockoff', 'AnnoKn', 'GhostKnockoff','AnnoGK'),]

methods <- c("Knockoff", "AnnoKn", "GhostKnockoff",
             "AnnoGK", "GhostKnockoff M=5", "AnnoGK M=5")

colors <- c("Knockoff" = "#1b9e77",
            "AnnoKn" = "#d95f02",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#e7298a",
            "GhostKnockoff M=5" = "#66a61e",
            "AnnoGK M=5" = "#e6ab02")

shapes <- c("Knockoff" = 16,
            "AnnoKn" = 17,
            "GhostKnockoff" = 15,
            "AnnoGK" = 3,
            "GhostKnockoff M=5" = 8,
            "AnnoGK M=5" = 4)

p1 <- ggplot(plot_data, aes(x = alpha, y = power, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Target FDR", y = "Power") +
  scale_color_manual(values = colors, breaks = methods) +
  scale_shape_manual(values = shapes, breaks = methods) +
  theme_minimal()

p2 <- ggplot(plot_data, aes(x = alpha, y = fdr, color = method, shape = method)) +
  geom_line() +
  geom_point() +
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_segment(
    aes(x = 0.1, y = 0.1, xend = 0.4, yend = 0.4), 
    linetype = "solid", 
    color = "black", 
    linewidth = 1
  ) +
  labs(x = "Target FDR", y = "FDR") +
  scale_color_manual(values = colors, breaks = methods) +
  scale_shape_manual(values = shapes, breaks = methods) +
  theme_minimal() +
  coord_cartesian(xlim = c(0.1, 0.4), ylim = c(0, 0.4))

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)












#################################### continuous ####################################



p = 600


library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.4, 0.1, by = -0.05)
len = length(alphalist)
heri_list <- c(0.05, 0.1, 0.2)


all_plot_data <- data.frame()

for (h in heri_list) {
  
  power_list <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))
  fdr_list   <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))
  
  removelist <- c()
  
  for (i in 1:iteration) {
    path <- paste0(repro_path("annogk_results", "heri_final_10_"),
                   h, "_n_", 5000, "_p_", p, "_", i, ".RData")
    if (!file.exists(path)) {
      removelist <- c(removelist, i)
      next
    }
    load(path)
    
    for (j in 1:11) {
      power_list[[j]][i, ] <- result[[paste0("power", j)]]
      fdr_list[[j]][i, ]   <- result[[paste0("fdr", j)]]
    }
  }
  
  power_means <- unlist(lapply(power_list, function(x) colMeans(x, na.rm = TRUE)))
  fdr_means   <- unlist(lapply(fdr_list, function(x) colMeans(x, na.rm = TRUE)))
  
  name_method <- c('Knockoff', 
                   'AnnoKn-lite', 
                   'AnnoKn', 
                   'GhostKnockoff', 
                   'AnnoGK-lite',
                   'AnnoGK',
                   'AnnoGK-dss', 
                   'GhostKnockoff M=5', 
                   'AnnoGK-simple M=5',
                   'AnnoGK M=5', 
                   'AnnoGK-dss M=5')
  
  methods <- c("Knockoff", "AnnoKn", "GhostKnockoff",
               "AnnoGK")
  
  method <- rep(name_method, each = len)
  tmp_df <- data.frame(
    alpha = rep(alphalist, 11),
    method = method,
    power = power_means,
    fdr = fdr_means,
    heritability = paste0("h\u00B2 = ", h)   # facet label
  )
 
  tmp_df <- tmp_df[tmp_df$method %in% methods, ]
  
  tmp_df_long <- bind_rows(
    tmp_df %>% transmute(alpha, value = power, method, heritability, measure = "Power"),
    tmp_df %>% transmute(alpha, value = fdr,   method, heritability, measure = "FDR")
  )
  
  all_plot_data <- rbind(all_plot_data, tmp_df_long)
}

all_plot_data$measure <- factor(all_plot_data$measure, levels = c("Power", "FDR"))

colors <- c("Knockoff" = "#e6ab02",
            "AnnoKn" = "#F8766D",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#4DAF4A",
            "AnnoKn-lite" = "#619CFF",
            "AnnoGK-lite" = "#D55E00")


shapes <- c("Knockoff" = 16,
            "AnnoKn" = 8,
            "GhostKnockoff" = 15,
            "AnnoGK" = 3,
            "AnnoKn-lite" = 5,
            "AnnoGK-lite" = 4)

linetypes <- c("Knockoff" = "solid",
               "AnnoKn" = "dashed",
               "GhostKnockoff" = "23",
               "AnnoGK" = "21",
               "AnnoKn-lite" = "dotted",
               "AnnoGK-lite" = "24")

plot_data_with_limits <- all_plot_data %>%
  mutate(ymax = ifelse(measure == "FDR", 0.4, 1))

g <- ggplot(all_plot_data, aes(x = alpha, y = value, color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  facet_grid(measure ~ heritability, switch = "y", scales = "free_y") +  
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

fdr_data <- subset(all_plot_data, measure == "FDR")
g <- g + geom_line(
  data = fdr_data,
  aes(x = alpha, y = alpha),
  inherit.aes = FALSE,
  linetype = "solid",
  color = "black",
  linewidth = 1
)

g <- g + geom_blank(
  data = plot_data_with_limits,
  aes(y = ymax)
)

print(g)

ggsave(repro_path("plots", "AnnoGK_simu_strong_p600.pdf"),
       g, 
       width = 9, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)







plot(beta,
     type = "h",                    
     lwd = 2,                      
     col = ifelse(beta != 0, "red", "grey"),
     xlab = "Index",
     ylab = "True coefficient")

box(bty = "n")





#################################### binary ####################################



p = 300


library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.4, 0.1, by = -0.05)
len = length(alphalist)
heri_list <- c(0.05, 0.1, 0.2)


all_plot_data <- data.frame()

for (h in heri_list) {
  
  power_list <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))
  fdr_list   <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))
  
  removelist <- c()
  
  for (i in 1:iteration) {
    path <- paste0(repro_path("annogk_results", "heri_final_10_binary_"),
                   h, "_n_", 5000, "_p_", p, "_", i, ".RData")
    if (!file.exists(path)) {
      removelist <- c(removelist, i)
      next
    }
    load(path)
    
    for (j in 1:11) {
      power_list[[j]][i, ] <- result[[paste0("power", j)]]
      fdr_list[[j]][i, ]   <- result[[paste0("fdr", j)]]
    }
  }
  
  power_means <- unlist(lapply(power_list, function(x) colMeans(x, na.rm = TRUE)))
  fdr_means   <- unlist(lapply(fdr_list, function(x) colMeans(x, na.rm = TRUE)))
  
  name_method <- c('Knockoff', 
                   'AnnoKn-lite', 
                   'AnnoKn', 
                   'GhostKnockoff', 
                   'AnnoGK-lite',
                   'AnnoGK',
                   'AnnoGK-dss', 
                   'GhostKnockoff M=5', 
                   'AnnoGK-simple M=5',
                   'AnnoGK M=5', 
                   'AnnoGK-dss M=5')
  
  methods <- c("Knockoff", "AnnoKn", "GhostKnockoff",
               "AnnoGK")
  
  method <- rep(name_method, each = len)
  tmp_df <- data.frame(
    alpha = rep(alphalist, 11),
    method = method,
    power = power_means,
    fdr = fdr_means,
    heritability = paste0("h\u00B2 = ", h)   # facet label
  )
  
  tmp_df <- tmp_df[tmp_df$method %in% methods, ]
  
  tmp_df_long <- bind_rows(
    tmp_df %>% transmute(alpha, value = power, method, heritability, measure = "Power"),
    tmp_df %>% transmute(alpha, value = fdr,   method, heritability, measure = "FDR")
  )
  
  all_plot_data <- rbind(all_plot_data, tmp_df_long)
}

all_plot_data$measure <- factor(all_plot_data$measure, levels = c("Power", "FDR"))

colors <- c("Knockoff" = "#e6ab02",
            "AnnoKn" = "#F8766D",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#4DAF4A",
            "AnnoKn-lite" = "#619CFF",
            "AnnoGK-lite" = "#D55E00")


# colors <- c('Adapt' = '#1b9e77', 
#             'Knockoff' = '#e6ab02', 
#             'AdaKn (GLM)' = '#d95f02', 
#             'AdaKn (GAM)' = '#00BFC4', 
#             'AdaKn (RF)' = '#e7298a',
#             'AdaKn (EM)' = '#C77CFF',
#             'AnnoKn' = '#F8766D', 
#             'AnnoKn-lite' = '#619CFF')

# shapes <- c('Adapt' = 11,
#             "Knockoff" = 16,
#             'AdaKn (GLM)' = 12,
#             'AdaKn (GAM)' = 13,
#             'AdaKn (RF)' = 17,
#             'AdaKn (EM)' = 7,
#             "AnnoKn" = 8,
#             "AnnoKn-lite" = 5)

shapes <- c("Knockoff" = 16,
            "AnnoKn" = 8,
            "GhostKnockoff" = 15,
            "AnnoGK" = 3,
            "AnnoKn-lite" = 5,
            "AnnoGK-lite" = 4)

# linetypes <- c('Adapt' = "82",
#                "Knockoff" = "solid",
#                'AdaKn (GLM)' = "41",
#                'AdaKn (GAM)' = "14",
#                "AdaKn (RF)" = "22",
#                'AdaKn (EM)' = "81",
#                "AnnoKn" = "dashed",
#                "AnnoKn-lite" = "dotted")

linetypes <- c("Knockoff" = "solid",
               "AnnoKn" = "dashed",
               "GhostKnockoff" = "23",
               "AnnoGK" = "21",
               "AnnoKn-lite" = "dotted",
               "AnnoGK-lite" = "24")

plot_data_with_limits <- all_plot_data %>%
  mutate(ymax = ifelse(measure == "FDR", 0.4, 1))

g <- ggplot(all_plot_data, aes(x = alpha, y = value, color = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  facet_grid(measure ~ heritability, switch = "y", scales = "free_y") +  
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

fdr_data <- subset(all_plot_data, measure == "FDR")
g <- g + geom_line(
  data = fdr_data,
  aes(x = alpha, y = alpha),
  inherit.aes = FALSE,
  linetype = "solid",
  color = "black",
  linewidth = 1
)

g <- g + geom_blank(
  data = plot_data_with_limits,
  aes(y = ymax)
)

print(g)

ggsave(repro_path("plots", "AnnoGK_simu_strong_binary_p300.pdf"),
       g, 
       width = 9, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)






######################### running time


p = 1000


library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.4, 0.1, by = -0.05)
len = length(alphalist)


all_plot_data <- data.frame()

h = 0.1

removelist <- c()
time1 <- time2 <- time3 <- time4 <- time5 <- time6 <- numeric()
time7 <- time8 <- time9 <- time10 <- time11 <- numeric()

for (i in 1:iteration) {
  path <- paste0(repro_path("annogk_results", "heri_final_10_"),
                 h, "_n_", 5000, "_p_", p, "_", i, ".RData")
  if (!file.exists(path)) {
    removelist <- c(removelist, i)
    next
  }
  load(path)
  
  time1 <- append(time1, result$time1)
  time2 <- append(time2, result$time2)
  time3 <- append(time3, result$time3)
  time4 <- append(time4, result$time4)
  time5 <- append(time5, result$time5)
  time6 <- append(time6, result$time6)
  time7 <- append(time7, result$time7)
  time8 <- append(time8, result$time8)
  time9 <- append(time9, result$time9)
  time10 <- append(time10, result$time10)
  time11 <- append(time11, result$time11)
}


name_method <- c('Knockoff', 
                 'AnnoKn-lite', 
                 'AnnoKn', 
                 'GhostKnockoff', 
                 'AnnoGK-lite',
                 'AnnoGK',
                 'AnnoGK-dss', 
                 'GhostKnockoff M=5', 
                 'AnnoGK-simple M=5',
                 'AnnoGK M=5', 
                 'AnnoGK-dss M=5')

methods <- c("Knockoff", "AnnoKn", "GhostKnockoff",
             "AnnoGK")

df <- data.frame(
  time = c(time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11),
  method = factor(rep(name_method, each = 100), levels = name_method)
)

tmp_df <- df[df$method %in% methods, ]


colors <- c("Knockoff" = "#e6ab02",
            "AnnoKn" = "#F8766D",
            "GhostKnockoff" = "#7570b3",
            "AnnoGK" = "#4DAF4A",
            "AnnoKn-lite" = "#619CFF",
            "AnnoGK-lite" = "#D55E00")

g <- ggplot(tmp_df, aes(x = method, y = time, fill = method)) +
  geom_boxplot(outlier.size = 1.2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    name = "Time (seconds)", 
    labels = function(x) paste0(x, " s")
  ) +
  labs(x = "Method") +
  theme_bw(base_size = 14) +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

g


ggsave(repro_path("plots", "AnnoGK_simu_time_h0.1_p1000.pdf"),
       g, 
       width = 8, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)









