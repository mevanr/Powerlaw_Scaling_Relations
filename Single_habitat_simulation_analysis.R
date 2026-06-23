# ============================================================================
# HABITAT DATA VISUALIZATION FOR POWER-LAW SIMULATIONS - R VERSION
# 6-HABITAT AND 12-HABITAT CASES DISPLAYED ON R SCREEN
# ============================================================================

rm(list = ls())

packages <- c(
  "tidyverse", "ggplot2", "gridExtra", "cowplot", "viridis",
  "MASS", "lme4", "ggpubr", "patchwork", "RColorBrewer"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# GLOBAL PARAMETERS
# ============================================================================

n_sim <- 1
n_per_habitat <- 36
multi_habitat <- TRUE

beta_values <- c(0.6, 1.0, 1.4)
alpha_true <- 1.5

sd_beta <- 0.40
sd_alpha <- 0.40
cor_alpha_beta <- -0.3

x_min <- 0.5
x_max <- 4
x_sd_habitat <- 1.2

sigma_x_me <- 0.20
sigma_y_me <- 0.40
sigma_proc <- 0.30
heteroscedastic <- FALSE

out_dir <- "E:/LMEM_Seagate/Habitat_mix_plots"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("\n📁 Output directory:", out_dir, "\n")

# ============================================================================
# THEME WITHOUT GRID LINES
# ============================================================================

theme_nogrid <- function() {
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

generate_correlated_random_effects <- function(n, sd_a, sd_b, corr) {
  Sigma <- matrix(c(sd_a^2, corr * sd_a * sd_b,
                    corr * sd_a * sd_b, sd_b^2), 2, 2)
  MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
}

ols_slope <- function(log_x, log_y) {
  fit <- lm(log_y ~ log_x)
  unname(coef(fit)[2])
}

theil_sen_slope <- function(log_x, log_y) {
  n <- length(log_x)
  slopes <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (abs(log_x[j] - log_x[i]) > 1e-10) {
        slopes <- c(slopes, (log_y[j] - log_y[i]) / (log_x[j] - log_x[i]))
      }
    }
  }
  if (length(slopes) == 0) return(NA)
  median(slopes)
}

# ============================================================================
# DATA GENERATION FUNCTION
# ============================================================================

generate_habitat_data <- function(beta_true, sim_id) {

  set.seed(42 + sim_id * 1000 + beta_true * 100 + n_habitats)

  cat(sprintf("\n  Generating data for H=%d, β=%.1f, simulation %d",
              n_habitats, beta_true, sim_id))

  habitat <- rep(1:n_habitats, each = n_per_habitat)
  N <- length(habitat)

  random_effects <- generate_correlated_random_effects(
    n_habitats, sd_alpha, sd_beta, cor_alpha_beta
  )

  beta_h <- beta_true + random_effects[, 2]
  alpha_h <- pmax(alpha_true + random_effects[, 1], 0.1)

  log_x_true <- numeric(N)

  for (h in 1:n_habitats) {
    idx <- which(habitat == h)
    log_mean <- log(x_means_habitat[h])
    log_sd <- log(x_sd_habitat + 1) / 2

    log_x_true[idx] <- rnorm(
      length(idx),
      mean = log_mean,
      sd = log_sd
    )
  }

  log_y_true <- numeric(N)

  for (i in 1:N) {
    h <- habitat[i]
    log_y_true[i] <- log(alpha_h[h]) +
      beta_h[h] * log_x_true[i] +
      rnorm(1, 0, sigma_proc)
  }

  if (heteroscedastic) {
    x_weights <- exp(log_x_true) / max(exp(log_x_true))
    y_weights <- exp(log_y_true) / max(exp(log_y_true))

    log_x_obs <- log_x_true + rnorm(N, 0, sigma_x_me * x_weights)
    log_y_obs <- log_y_true + rnorm(N, 0, sigma_y_me * y_weights)

    se_x <- sigma_x_me * x_weights
    se_y <- sigma_y_me * y_weights
  } else {
    log_x_obs <- log_x_true + rnorm(N, 0, sigma_x_me)
    log_y_obs <- log_y_true + rnorm(N, 0, sigma_y_me)

    se_x <- rep(sigma_x_me, N)
    se_y <- rep(sigma_y_me, N)
  }

  df <- data.frame(
    sim_id = sim_id,
    beta_true = beta_true,
    n_habitats = n_habitats,
    habitat = factor(habitat),
    log_x_true = log_x_true,
    log_y_true = log_y_true,
    log_x_obs = log_x_obs,
    log_y_obs = log_y_obs,
    se_x = se_x,
    se_y = se_y,
    x_true = exp(log_x_true),
    y_true = exp(log_y_true),
    x_obs = exp(log_x_obs),
    y_obs = exp(log_y_obs),
    alpha_true_habitat = alpha_h[habitat],
    beta_true_habitat = beta_h[habitat]
  )

  df$beta_ols <- ols_slope(df$log_x_obs, df$log_y_obs)
  df$beta_theilsen <- theil_sen_slope(df$log_x_obs, df$log_y_obs)

  return(df)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

plot_single_habitat_data <- function(df, out_dir) {

  sim_id <- unique(df$sim_id)
  beta_true <- unique(df$beta_true)
  H <- unique(df$n_habitats)

  n_hab <- length(unique(df$habitat))
  colors <- viridis::viridis(n_hab)
  names(colors) <- levels(df$habitat)

  p1 <- ggplot(df, aes(x = log_x_obs, y = log_y_obs, color = habitat)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colors) +
    labs(
      title = paste("Simulation", sim_id, "- True β =", beta_true, "- H =", H),
      subtitle = "Habitat-specific data with true slopes (solid) and fitted lines (dashed)",
      x = "log(X)",
      y = "log(Y)"
    ) +
    theme_nogrid() +
    theme(legend.position = "bottom")

  for (h in levels(df$habitat)) {
    hab_data <- df[df$habitat == h, ]

    if (nrow(hab_data) > 0) {
      alpha_h <- log(unique(hab_data$alpha_true_habitat))
      beta_h <- unique(hab_data$beta_true_habitat)

      p1 <- p1 +
        geom_abline(
          intercept = alpha_h,
          slope = beta_h,
          color = colors[h],
          linewidth = 1.2,
          alpha = 0.8
        )

      fit <- lm(log_y_obs ~ log_x_obs, data = hab_data)

      p1 <- p1 +
        geom_abline(
          intercept = coef(fit)[1],
          slope = coef(fit)[2],
          color = colors[h],
          linetype = "dashed",
          alpha = 0.5
        )
    }
  }

  p2 <- ggplot(df, aes(x = habitat, y = log_x_obs, fill = habitat)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = colors) +
    labs(
      title = "Distribution of log(X) by Habitat",
      x = "Habitat",
      y = "log(X)"
    ) +
    theme_nogrid() +
    theme(legend.position = "none")

  methods_df <- data.frame(
    method = c("OLS", "Theil-Sen", paste("Habitat", 1:min(3, n_hab), "True")),
    estimate = c(
      unique(df$beta_ols),
      unique(df$beta_theilsen),
      sapply(1:min(3, n_hab), function(i) {
        unique(df$beta_true_habitat[df$habitat == i])[1]
      })
    )
  )

  p3 <- ggplot(methods_df, aes(x = method, y = estimate, fill = method)) +
    geom_col(alpha = 0.7) +
    geom_hline(
      yintercept = beta_true,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    geom_text(aes(label = round(estimate, 3)), vjust = -0.5) +
    labs(
      title = "Method Comparison",
      x = "",
      y = "β Estimate"
    ) +
    theme_nogrid() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  fit <- lm(log_y_obs ~ log_x_obs, data = df)
  df$residuals <- resid(fit)

  p4 <- ggplot(df, aes(x = log_x_obs, y = residuals, color = habitat)) +
    geom_point(alpha = 0.6) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red"
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = paste("Residuals (R² =", round(summary(fit)$r.squared, 3), ")"),
      x = "log(X)",
      y = "Residuals"
    ) +
    theme_nogrid() +
    theme(legend.position = "none")

  p5 <- ggplot(df, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(
      title = "Q-Q Plot of Residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_nogrid()

  summary_text <- paste(
    "SUMMARY STATISTICS:\n",
    sprintf("Sample Size: %d\n", nrow(df)),
    sprintf("Habitats: %d\n", n_hab),
    sprintf("σx = %.2f, σy = %.2f, σproc = %.2f\n",
            sigma_x_me, sigma_y_me, sigma_proc),
    "Habitat X Means:\n",
    paste(sprintf("  H%d: %.1f", 1:n_habitats, x_means_habitat), collapse = "\n"),
    "\nRandom Effects:\n",
    sprintf("  SD(α) = %.3f\n", sd(df$alpha_true_habitat)),
    sprintf("  SD(β) = %.3f\n", sd(df$beta_true_habitat)),
    sprintf("  Cor(α,β) = %.3f", cor(df$alpha_true_habitat, df$beta_true_habitat))
  )

  top_row <- plot_grid(p1, p2, ncol = 2, rel_widths = c(2, 1))
  middle_row <- plot_grid(p3, ncol = 1)
  bottom_row <- plot_grid(p4, p5, ncol = 2)

  final_plot <- plot_grid(
    top_row,
    middle_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(1.2, 0.8, 1)
  )

  title <- ggdraw() +
    draw_label(
      paste("Habitat-Specific Data Analysis - H =", H, "- Simulation", sim_id),
      fontface = "bold",
      size = 16
    )

  final_plot <- plot_grid(
    title,
    final_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  print(final_plot)

  filename <- file.path(
    out_dir,
    sprintf("habitat_H%d_sim%d_beta%.2f.png", H, sim_id, beta_true)
  )

  tryCatch({
    ggsave(filename, final_plot, width = 16, height = 14, dpi = 150, units = "in")

    if (file.exists(filename)) {
      cat(sprintf("\n    ✓ Plot saved: %s (%d bytes)", filename, file.size(filename)))
    } else {
      cat(sprintf("\n    ❌ Plot not saved: %s", filename))
    }
  }, error = function(e) {
    cat(sprintf("\n    ❌ Error saving plot: %s", e$message))
  })

  txt_file <- file.path(
    out_dir,
    sprintf("summary_H%d_sim%d_beta%.2f.txt", H, sim_id, beta_true)
  )

  writeLines(summary_text, txt_file)

  if (file.exists(txt_file)) {
    cat(sprintf("\n    ✓ Text summary saved: %s", txt_file))
  }
}

plot_habitat_comparison <- function(all_dfs, out_dir, H) {

  all_data <- bind_rows(all_dfs)

  p <- ggplot(all_data, aes(x = log_x_obs, y = log_y_obs, color = factor(beta_true))) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    facet_wrap(~beta_true, ncol = 3, labeller = label_both) +
    scale_color_viridis_d() +
    labs(
      title = paste("All Simulations Combined - H =", H),
      x = "log(X)",
      y = "log(Y)",
      color = "True β"
    ) +
    theme_nogrid() +
    theme(legend.position = "bottom")

  print(p)

  filename <- file.path(out_dir, sprintf("all_simulations_comparison_H%d.png", H))
  ggsave(filename, p, width = 15, height = 5, dpi = 150)

  cat("\n  ✓ Comparison plot saved:", filename)
}

plot_habitat_summary <- function(all_dfs, out_dir, H) {

  all_data <- bind_rows(all_dfs)

  p1 <- ggplot(all_data, aes(x = habitat, y = beta_true_habitat, fill = habitat)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(
      yintercept = beta_values,
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
    scale_fill_viridis_d() +
    labs(
      title = "Distribution of True Slopes by Habitat",
      x = "Habitat",
      y = "True β"
    ) +
    theme_nogrid() +
    theme(legend.position = "none")

  x_means <- all_data %>%
    group_by(habitat) %>%
    summarise(mean_x = mean(x_true), .groups = "drop")

  p2 <- ggplot(x_means, aes(x = habitat, y = mean_x, fill = habitat)) +
    geom_col(alpha = 0.7) +
    geom_hline(
      yintercept = mean(x_means$mean_x),
      linetype = "dashed",
      color = "red"
    ) +
    scale_fill_viridis_d() +
    labs(
      title = "Habitat-Specific X Means",
      x = "Habitat",
      y = "Mean X"
    ) +
    theme_nogrid() +
    theme(legend.position = "none")

  estimates <- all_data %>%
    group_by(sim_id, beta_true) %>%
    summarise(
      ols = first(beta_ols),
      theilsen = first(beta_theilsen),
      .groups = "drop"
    )

  p3 <- ggplot(estimates, aes(x = ols, y = theilsen, color = factor(beta_true))) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "red"
    ) +
    scale_color_viridis_d() +
    labs(
      title = "OLS vs Theil-Sen Estimates",
      x = "OLS Estimate",
      y = "Theil-Sen Estimate",
      color = "True β"
    ) +
    theme_nogrid() +
    theme(legend.position = "bottom")

  errors <- estimates %>%
    mutate(error = ols - beta_true)

  p4 <- ggplot(errors, aes(x = error, fill = factor(beta_true))) +
    geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    scale_fill_viridis_d() +
    labs(
      title = "Distribution of OLS Estimation Errors",
      x = "Error (OLS - True)",
      y = "Frequency",
      fill = "True β"
    ) +
    theme_nogrid() +
    theme(legend.position = "bottom")

  top_row <- plot_grid(p1, p2, ncol = 2)
  bottom_row <- plot_grid(p3, p4, ncol = 2)

  final_plot <- plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(1, 1.2)
  )

  title <- ggdraw() +
    draw_label(
      paste("Summary Statistics Across All Simulations - H =", H),
      fontface = "bold",
      size = 16
    )

  final_plot <- plot_grid(
    title,
    final_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  print(final_plot)

  filename <- file.path(out_dir, sprintf("summary_statistics_H%d.png", H))
  ggsave(filename, final_plot, width = 14, height = 12, dpi = 150)

  cat("\n  ✓ Summary plot saved:", filename)
}

create_data_report <- function(all_dfs, out_dir, H) {

  all_data <- bind_rows(all_dfs)

  report_file <- file.path(out_dir, sprintf("simulation_report_H%d.txt", H))
  sink(report_file)

  cat("========================================================\n")
  cat("HABITAT DATA SIMULATION REPORT\n")
  cat("========================================================\n\n")

  cat("SIMULATION PARAMETERS:\n")
  cat("---------------------\n")
  cat("Number of simulations per β:", n_sim, "\n")
  cat("Observations per habitat:", n_per_habitat, "\n")
  cat("Number of habitats:", H, "\n")
  cat("True β values:", paste(beta_values, collapse = ", "), "\n")
  cat("True α:", alpha_true, "\n\n")

  cat("Random Effects:\n")
  cat("  SD(β):", sd_beta, "\n")
  cat("  SD(α):", sd_alpha, "\n")
  cat("  Cor(α,β):", cor_alpha_beta, "\n\n")

  cat("Error Parameters:\n")
  cat("  σx:", sigma_x_me, "\n")
  cat("  σy:", sigma_y_me, "\n")
  cat("  σproc:", sigma_proc, "\n")
  cat("  Heteroscedastic:", heteroscedastic, "\n\n")

  cat("Habitat X Means:\n")
  for (i in 1:H) {
    cat(sprintf("  Habitat %d: %.1f\n", i, x_means_habitat[i]))
  }

  cat("\nSUMMARY STATISTICS:\n")
  cat("------------------\n")
  cat("Total observations:", nrow(all_data), "\n")
  cat("Total simulations:", length(unique(all_data$sim_id)), "\n\n")

  for (beta_val in beta_values) {
    beta_data <- all_data[all_data$beta_true == beta_val, ]

    cat(sprintf("\nTrue β = %.1f:\n", beta_val))
    cat(sprintf("  Number of simulations: %d\n", length(unique(beta_data$sim_id))))
    cat(sprintf("  Mean OLS estimate: %.3f\n", mean(beta_data$beta_ols, na.rm = TRUE)))
    cat(sprintf("  SD of OLS estimates: %.3f\n", sd(beta_data$beta_ols, na.rm = TRUE)))
    cat(sprintf("  Mean Theil-Sen estimate: %.3f\n", mean(beta_data$beta_theilsen, na.rm = TRUE)))
    cat(sprintf("  SD of Theil-Sen estimates: %.3f\n", sd(beta_data$beta_theilsen, na.rm = TRUE)))
  }

  sink()

  csv_file <- file.path(out_dir, sprintf("all_simulation_data_H%d.csv", H))
  write.csv(all_data, csv_file, row.names = FALSE)

  cat("\n  ✓ Report saved:", report_file)
  cat("\n  ✓ CSV data saved:", csv_file)
}

# ============================================================================
# MAIN EXECUTION: RUN 6-HABITAT AND 12-HABITAT CASES
# ============================================================================

cat("\n")
cat("========================================================\n")
cat("HABITAT DATA SIMULATION AND VISUALIZATION - R VERSION\n")
cat("6-HABITAT AND 12-HABITAT CASES\n")
cat("========================================================\n")

# ============================================================================
# ADDED: OBJECTS TO SAVE ALL RESULTS ACROSS BOTH HABITAT CASES
# ============================================================================

overall_all_dfs <- list()

for (H_case in c(6, 12)) {

  n_habitats <- H_case

  if (n_habitats == 6) {
    x_means_habitat <- c(1.0, 2.5, 4.0, 5.5, 7.0, 8.5)
  }

  if (n_habitats == 12) {
    x_means_habitat <- c(
      1.0, 1.8, 2.6, 3.4, 4.2, 5.0,
      5.8, 6.6, 7.4, 8.2, 9.0, 9.8
    )
  }

  cat("\n")
  cat("========================================================\n")
  cat(sprintf("RUNNING %d-HABITAT CASE\n", n_habitats))
  cat("========================================================\n")

  cat("\n📊 Simulation Parameters:\n")
  cat(sprintf("  • Simulations per β: %d\n", n_sim))
  cat(sprintf("  • Habitats: %d with %d observations each\n", n_habitats, n_per_habitat))
  cat(sprintf("  • True β values: %s\n", paste(beta_values, collapse = ", ")))
  cat(sprintf("  • Random slopes SD: %.2f\n", sd_beta))
  cat(sprintf("  • Random intercepts SD: %.2f\n", sd_alpha))
  cat(sprintf("  • Correlation (α-β): %.2f\n", cor_alpha_beta))
  cat(sprintf("  • Habitat X means: %s\n", paste(x_means_habitat, collapse = ", ")))
  cat(sprintf("  • Output directory: %s\n", out_dir))
  cat(sprintf("  • Working directory: %s\n", getwd()))

  all_dfs <- list()

  for (beta_true in beta_values) {
    cat(sprintf("\n▶ Processing β = %.1f", beta_true))

    for (sim_id in 1:n_sim) {
      df <- generate_habitat_data(beta_true, sim_id)
      all_dfs[[length(all_dfs) + 1]] <- df

      sim_file <- file.path(
        out_dir,
        sprintf("data_H%d_sim%d_beta%.2f.csv", n_habitats, sim_id, beta_true)
      )

      write.csv(df, sim_file, row.names = FALSE)

      if (file.exists(sim_file)) {
        cat(sprintf("\n    ✓ Data saved: %s", sim_file))
      }

      plot_single_habitat_data(df, out_dir)
    }
  }

  cat(sprintf("\n\n✅ Generated %d simulations total for H=%d\n",
              length(all_dfs), n_habitats))

  cat("\n\n📊 Creating comparison plots...\n")
  plot_habitat_comparison(all_dfs, out_dir, n_habitats)
  plot_habitat_summary(all_dfs, out_dir, n_habitats)

  cat("\n\n📝 Creating data report...\n")
  create_data_report(all_dfs, out_dir, n_habitats)

  # ==========================================================================
  # ADDED: STORE THIS HABITAT CASE FOR FINAL COMBINED SAVING
  # ==========================================================================
  overall_all_dfs[[as.character(n_habitats)]] <- bind_rows(all_dfs)
}

# ============================================================================
# ADDED: SAVE COMBINED RESULTS ACROSS BOTH 6-HABITAT AND 12-HABITAT CASES
# ============================================================================

overall_data <- bind_rows(overall_all_dfs, .id = "H_case_saved")

overall_csv <- file.path(out_dir, "all_simulation_data_H6_and_H12_combined.csv")
write.csv(overall_data, overall_csv, row.names = FALSE)

overall_rds <- file.path(out_dir, "all_simulation_data_H6_and_H12_combined.rds")
saveRDS(overall_data, overall_rds)

workspace_file <- file.path(out_dir, "habitat_simulation_workspace.RData")
save.image(workspace_file)

cat("\n  ✓ Combined CSV data saved:", overall_csv)
cat("\n  ✓ Combined RDS data saved:", overall_rds)
cat("\n  ✓ Full R workspace saved:", workspace_file)

cat("\n\n")
cat("========================================================\n")
cat("✅ ALL 6-HABITAT AND 12-HABITAT SIMULATIONS COMPLETE\n")
cat("========================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("========================================================\n")