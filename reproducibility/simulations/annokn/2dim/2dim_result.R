# Paths are configured by reproducibility/run.R.
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)

cairo_pdf(repro_path("revision_plots", "AnnoKn_2dim_type.pdf"),
          width = 10, height = 10, bg = "white")

par(mfrow = c(2, 2))

type_all = c(1, 0.5, 2, 'in')

for(i in 1:4){
  type = type_all[i]
  
  
  nonzero = c()
  
  cor_hypothesis <- expand.grid(x = 1:40, y = 1:40)
  cor_hypothesis$index <- 1:nrow(cor_hypothesis)
  # causal_region <- subset(cor_hypothesis, x <= 20 & y <= 20)
  
  if(type == 1){
    causal_region <- subset(cor_hypothesis, x + y <= 30)
  }else if(type == 2){
    causal_region <- subset(cor_hypothesis, x^2 + y^2 <= 600)
  }else if(type == 0.5){
    causal_region <- subset(cor_hypothesis, sqrt(x) + sqrt(y) <= sqrt(50))
  }else if(type == 'in'){
    causal_region <- subset(cor_hypothesis, x <= 20 & y <= 20)
  }
  
  causal_region$prob <- 1 / (causal_region$x + causal_region$y)
  #causal_region$prob <- 1/((causal_region$x)^d + (causal_region$y)^d)
  causal_region$prob <- causal_region$prob / sum(causal_region$prob)  
  selected_rows <- sample(1:nrow(causal_region), size = k, prob = causal_region$prob)
  nonzero <- causal_region$index[selected_rows[order(selected_rows)]]
  
  
  ############################ making the plot
  
  indi <- numeric(p)
  indi[nonzero] = 1
  col_vec <- ifelse(indi == 1, "black", "grey")
  
  plot(x = cor_hypothesis$x,
       y = cor_hypothesis$y,
       col = col_vec,
       pch = 16,
       xlab = "r",
       ylab = "s",
       cex.lab = 1.5,
       main = paste("Scenario", i))
}

dev.off()









####################### result ######################


library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.3,0.05,-0.05)
amp = 25
binary = ''
cortype = 3
type_all = c(1, 0.5, 2, 'in')

len = length(alphalist)

all_plot_data <- data.frame()

for (type_i in 1:4) {
  type = type_all[type_i]
  
  power_list <- lapply(1:7, function(x) matrix(0, nrow = iteration, ncol = len))
  fdr_list   <- lapply(1:7, function(x) matrix(0, nrow = iteration, ncol = len))
  
  removelist <- c()
  
  for (i in 1:iteration) {
    path = paste0(repro_path("annokn_2d_results", "result_"),
                  binary, "_amp_",amp, "_type_", type, "_cor_", cortype, "_", i,".RData")
    
    
    if (!file.exists(path)) {
      removelist <- c(removelist, i)
      next
    }
    load(path)
    
    for (j in 1:7) {
      power_list[[j]][i, ] <- result[[paste0("power", j)]]
      fdr_list[[j]][i, ]   <- result[[paste0("fdr", j)]]
    }
  }
  
  power_means <- unlist(lapply(power_list, function(x) colMeans(x, na.rm = TRUE)))
  fdr_means   <- unlist(lapply(fdr_list, function(x) colMeans(x, na.rm = TRUE)))
  
  name_method <- c('Knockoff', 
                   'AdaKn (GLM)', 
                   'AdaKn (GAM)', 
                   'AdaKn (RF)', 
                   'AdaKn (EM)',
                   'AnnoKn',
                   'AnnoKn-lite')
  
  methods <- c('Knockoff', 
               # 'AdaKn (GLM)', 
               # 'AdaKn (GAM)', 
               'AdaKn (RF)', 
               'AdaKn (EM)',
               'AnnoKn',
               'AnnoKn-lite')
  
  method <- rep(name_method, each = len)
  tmp_df <- data.frame(
    alpha = rep(alphalist, 7),
    method = method,
    power = power_means,
    fdr = fdr_means,
    type = paste("Scenario", type_i)   # facet label
  )
  
  tmp_df <- tmp_df[tmp_df$method %in% methods, ]
  
  tmp_df_long <- bind_rows(
    tmp_df %>% transmute(alpha, value = power, method, type, measure = "Power"),
    tmp_df %>% transmute(alpha, value = fdr,   method, type, measure = "FDR")
  )
  
  all_plot_data <- rbind(all_plot_data, tmp_df_long)
}

all_plot_data$measure <- factor(all_plot_data$measure, levels = c("Power", "FDR"))




colors <- c('Adapt' = '#1b9e77', 
            'Knockoff' = '#e6ab02', 
            'AdaKn (GLM)' = '#d95f02', 
            'AdaKn (GAM)' = '#00BFC4', 
            'AdaKn (RF)' = '#e7298a',
            'AdaKn (EM)' = '#C77CFF',
            'AnnoKn' = '#F8766D', 
            'AnnoKn-lite' = '#619CFF')

shapes <- c('Adapt' = 11,
            "Knockoff" = 16,
            'AdaKn (GLM)' = 12,
            'AdaKn (GAM)' = 13,
            'AdaKn (RF)' = 17,
            'AdaKn (EM)' = 7,
            "AnnoKn" = 8,
            "AnnoKn-lite" = 5)

linetypes <- c('Adapt' = "82",
               "Knockoff" = "solid",
               'AdaKn (GLM)' = "41",
               'AdaKn (GAM)' = "14",
               "AdaKn (RF)" = "22",
               'AdaKn (EM)' = "81",
               "AnnoKn" = "dashed",
               "AnnoKn-lite" = "dotted")




plot_data_with_limits <- all_plot_data %>%
  mutate(ymax = ifelse(measure == "FDR", 0.3, 0.8))

g <- ggplot(all_plot_data, aes(x = alpha, y = value, color = method, 
                               shape = method, linetype = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  facet_grid(measure ~ type, switch = "y", scales = "free_y") +  
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


plot_path = paste0(repro_path("revision_plots", "AnnoKn_2dim_"), binary, "_", cortype, ".pdf")

ggsave(plot_path, 
       g, 
       width = 10, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)




################## running time:





library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.3,0.05,-0.05)
amp = 25
binary = ''
type = 1


time1 <- time2 <- time3 <- time4 <- time5 <- time6 <- time7 <- numeric()


removelist <- c()

for (i in 1:iteration) {
  path = paste0(repro_path("annokn_2d_timing_results", "result_"),
                binary, "_amp_",amp, "_type_", type, "_", i,".RData")
  
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
}


name_method <- c('Knockoff', 
                 'AdaKn (GLM)', 
                 'AdaKn (GAM)', 
                 'AdaKn (RF)', 
                 'AdaKn (EM)',
                 'AnnoKn',
                 'AnnoKn-lite')

methods <- c('Knockoff', 
             # 'AdaKn (GLM)', 
             # 'AdaKn (GAM)', 
             'AdaKn (RF)', 
             'AdaKn (EM)',
             'AnnoKn',
             'AnnoKn-lite')

df <- data.frame(
  time = c(time1, time2, time3, time4, time5, time6, time7),
  method = factor(rep(name_method, each = 100), levels = name_method)
)

tmp_df <- df[df$method %in% methods, ]


colors <- c('Adapt' = '#1b9e77', 
            'Knockoff' = '#e6ab02', 
            'AdaKn (GLM)' = '#d95f02', 
            'AdaKn (GAM)' = '#00BFC4', 
            'AdaKn (RF)' = '#e7298a',
            'AdaKn (EM)' = '#C77CFF',
            'AnnoKn' = '#F8766D', 
            'AnnoKn-lite' = '#619CFF')


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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

g

ggsave(repro_path("revision_plots", "AnnoKn_2dim_time_continuous.pdf"),
       g, 
       width = 6, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)




