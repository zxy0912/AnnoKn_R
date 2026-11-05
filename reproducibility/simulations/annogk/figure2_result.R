# Main manuscript Figure 2: standalone reader derived from result.R's continuous block.
# Statistical aggregation and method slots are retained; all 300 inputs are required.
# Run: Rscript reproducibility/run.R simulations/annogk/figure2_result.R
if (!exists("repro_path", mode = "function")) stop("Run this script through reproducibility/run.R; see reproducibility/README.md.", call. = FALSE)

p = 1000


library(ggplot2)
library(dplyr)

iteration = 100
alphalist <- seq(0.4, 0.1, by = -0.05)
len = length(alphalist)
heri_list <- c(0.05, 0.1, 0.2)


# Require all three settings and all 100 replicates before averaging.
input_files <- do.call(rbind, lapply(heri_list, function(h) {
  data.frame(h2 = h, replicate = seq_len(iteration),
             filename = paste0("heri_final_10_", h, "_n_5000_p_", p, "_",
                               seq_len(iteration), ".RData"))
}))
input_paths <- file.path(repro_path("annogk_results"), input_files$filename)
missing_files <- input_files$filename[!file.exists(input_paths)]
if (length(missing_files)) {
  stop("Figure 2 requires all 300 result files (100 per heritability). Missing ",
       length(missing_files), ": ", paste(head(missing_files, 10), collapse = ", "),
       if (length(missing_files) > 10) " ..." else "",
       ". Run continuous.R for the missing replicate IDs or configure annogk_results.",
       call. = FALSE)
}

all_plot_data <- data.frame()

for (h in heri_list) {

  power_list <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))
  fdr_list   <- lapply(1:11, function(x) matrix(0, nrow = iteration, ncol = len))


  for (i in 1:iteration) {
    path <- paste0(repro_path("annogk_results", "heri_final_10_"),
                   h, "_n_", 5000, "_p_", p, "_", i, ".RData")
    # Isolate loaded objects so saved workspaces cannot replace script settings.
    loaded <- new.env(parent = emptyenv())
    tryCatch(load(path, envir = loaded), error = function(e)
      stop("Cannot read Figure 2 result ", path, ": ", conditionMessage(e), call. = FALSE))
    if (!exists("result", envir = loaded, inherits = FALSE) || !is.list(loaded$result))
      stop("Expected a result list in ", path, call. = FALSE)
    result <- loaded$result
    fields <- c(paste0("power", 1:11), paste0("fdr", 1:11))
    for (field in fields) {
      value <- result[[field]]
      if (!is.numeric(value) || !is.null(dim(value)) || length(value) != len ||
          any(!is.finite(value)) || any(value < 0 | value > 1))
        stop("Invalid ", field, " in ", path,
             ": expected seven finite values between 0 and 1.", call. = FALSE)
    }

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
    text = element_text(family = "sans"),
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

# Save directly: no interactive/default device or Arial installation is required.
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(repro_path("plots", "AnnoGK_simu_strong_p1000.pdf"), g,
       width = 9, height = 5, units = "in", bg = "white", device = pdf_device)
write.table(all_plot_data,
            repro_path("plots", "AnnoGK_simu_strong_p1000_values.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
input_files$md5 <- unname(tools::md5sum(input_paths))
write.table(input_files,
            repro_path("plots", "AnnoGK_simu_strong_p1000_inputs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Figure 2 saved: continuous annotations, p=1000, 100 replicates per heritability.")
