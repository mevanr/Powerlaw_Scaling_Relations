# ============================================================================
# SINGLE-HABITAT POWER-LAW BENCHMARK
# ============================================================================
# Folder structure created exactly as:
#   habitat_visualization_single_output_n36_mex20_mey40
#   habitat_visualization_single_output_n36_mex30_mey30
#   habitat_visualization_single_output_n36_mex40_mey20
#   habitat_visualization_single_output_n72_mex20_mey40
#   habitat_visualization_single_output_n72_mex30_mey30
#   habitat_visualization_single_output_n72_mex40_mey20
#
# Each folder contains all three beta values: beta = 0.6, 1.0, 1.4.
# Results are written after each simulation replicate so progress is preserved.
# ============================================================================

rm(list = ls())

# ----------------------------------------------------------------------------
# 0. USER SETTINGS
# ----------------------------------------------------------------------------

QUICK_TEST <- FALSE
N_SIM_FINAL <- 500
N_SIM_TEST  <- 5
n_sim <- if (QUICK_TEST) N_SIM_TEST else N_SIM_FINAL

set.seed(20260623)

base_out_dir <- "D:/LMEM"
if (!dir.exists(base_out_dir)) dir.create(base_out_dir, recursive = TRUE)

# Manuscript single-habitat grid
n_values <- c(36, 72)
beta_values <- c(0.6, 1.0, 1.4)
error_configs <- tibble::tribble(
  ~sigma_x_me, ~sigma_y_me,
  0.20,        0.40,
  0.30,        0.30,
  0.40,        0.20
)

# Single-habitat latent power-law parameters
alpha_true <- 1.5
sigma_proc <- 0.30
x_min <- 0.5
x_max <- 8.5

# ----------------------------------------------------------------------------
# 1. PACKAGES
# ----------------------------------------------------------------------------

required_pkgs <- c(
  "tidyverse", "MASS", "quantreg", "glmnet", "mblm", "mcr", "lmodel2",
  "smatr", "e1071", "sandwich", "lmtest", "ggplot2", "readr"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it first with install.packages('%s')", pkg, pkg))
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(quantreg)
  library(glmnet)
  library(mblm)
  library(mcr)
  library(lmodel2)
  library(smatr)
  library(e1071)
  library(sandwich)
  library(lmtest)
})

# ----------------------------------------------------------------------------
# 2. EXACT SINGLE-HABITAT METHOD SET FROM THE MANUSCRIPT FIGURES
# ----------------------------------------------------------------------------
# Do not add or remove rows unless the manuscript method list changes.

METHOD_ORDER <- c(
  "ODR",
  "SMA (Major Axis)",
  "SMA",
  "RMA",
  "Theil-Sen",
  "Passing-Bablok",
  "Median",
  "Clustered SE",
  "OLS",
  "LMM (Random Slope)",
  "LMM (Random Intercept)",
  "Huber",
  "Trimmed Mean",
  "Quantile",
  "Ensemble",
  "Deming",
  "Elastic Net",
  "Lasso",
  "Ridge",
  "SVM Linear",
  "SVM RBF",
  "Adaptive Lasso"
)

METHOD_COLS <- c(
  beta_odr = "ODR",
  beta_sma_major_axis = "SMA (Major Axis)",
  beta_sma = "SMA",
  beta_rma = "RMA",
  beta_theilsen = "Theil-Sen",
  beta_passing_bablok = "Passing-Bablok",
  beta_median = "Median",
  beta_clustered_se = "Clustered SE",
  beta_ols = "OLS",
  beta_lmm_random_slope = "LMM (Random Slope)",
  beta_lmm_random_intercept = "LMM (Random Intercept)",
  beta_huber = "Huber",
  beta_trimmed_mean = "Trimmed Mean",
  beta_quantile = "Quantile",
  beta_ensemble = "Ensemble",
  beta_deming = "Deming",
  beta_elastic_net = "Elastic Net",
  beta_lasso = "Lasso",
  beta_ridge = "Ridge",
  beta_svm_linear = "SVM Linear",
  beta_svm_rbf = "SVM RBF",
  beta_adaptive_lasso = "Adaptive Lasso"
)

# ----------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ----------------------------------------------------------------------------

folder_name <- function(n_obs, sigma_x_me, sigma_y_me) {
  sprintf(
    "habitat_visualization_single_output_n%d_mex%02d_mey%02d",
    n_obs,
    round(100 * sigma_x_me),
    round(100 * sigma_y_me)
  )
}

scenario_label <- function(beta_true, n_obs, sigma_x_me, sigma_y_me) {
  sprintf(
    "beta=%s_n=%d_mex=%s_mey=%s",
    format(beta_true, trim = TRUE), n_obs,
    format(sigma_x_me, trim = TRUE), format(sigma_y_me, trim = TRUE)
  )
}

safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || is.na(x[1]) || is.nan(x[1]) || is.infinite(x[1])) return(NA_real_)
  as.numeric(x[1])
}

pairwise_slopes <- function(x, y) {
  n <- length(x)
  if (n < 2) return(numeric(0))
  idx <- utils::combn(seq_len(n), 2)
  dx <- x[idx[2, ]] - x[idx[1, ]]
  dy <- y[idx[2, ]] - y[idx[1, ]]
  slopes <- dy / dx
  slopes[is.finite(slopes) & abs(dx) > .Machine$double.eps]
}

median_slope <- function(x, y) {
  s <- pairwise_slopes(x, y)
  if (length(s) == 0) return(NA_real_)
  median(s, na.rm = TRUE)
}

trimmed_slope <- function(x, y, trim = 0.10) {
  s <- pairwise_slopes(x, y)
  if (length(s) == 0) return(NA_real_)
  mean(s, trim = trim, na.rm = TRUE)
}

# Unweighted major-axis/TLS slope
major_axis_slope <- function(x, y) {
  S <- stats::cov(cbind(x, y), use = "complete.obs")
  sxx <- S[1, 1]
  syy <- S[2, 2]
  sxy <- S[1, 2]
  if (!is.finite(sxy) || abs(sxy) < .Machine$double.eps) return(NA_real_)
  (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
}

# Weighted ODR/TLS after scaling by known measurement standard deviations.
# This is the single-habitat EIV correction used as the ODR row.
odr_known_error_slope <- function(x, y, sx, sy) {
  xs <- x / sx
  ys <- y / sy
  b_scaled <- major_axis_slope(xs, ys)
  if (!is.finite(b_scaled)) return(NA_real_)
  (sy / sx) * b_scaled
}

# Closed-form Deming with known error variance ratio lambda = sigma_y^2 / sigma_x^2.
deming_slope <- function(x, y, sx, sy) {
  lambda <- (sy^2) / (sx^2)
  S <- stats::cov(cbind(x, y), use = "complete.obs")
  sxx <- S[1, 1]
  syy <- S[2, 2]
  sxy <- S[1, 2]
  if (!is.finite(sxy) || abs(sxy) < .Machine$double.eps) return(NA_real_)
  (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 + 4 * lambda * sxy^2)) / (2 * sxy)
}

prediction_grid_slope <- function(model, predict_fun, x_min_obs, x_max_obs, grid_n = 100) {
  grid <- data.frame(log_x = seq(x_min_obs, x_max_obs, length.out = grid_n))
  pred <- tryCatch(predict_fun(model, grid), error = function(e) rep(NA_real_, nrow(grid)))
  if (all(is.na(pred))) return(NA_real_)
  safe_num(coef(lm(as.numeric(pred) ~ grid$log_x))[2])
}

fit_predict_component <- function(method, train_df, test_df) {
  tryCatch({
    if (method == "OLS") {
      fit <- lm(log_y ~ log_x, data = train_df)
      return(list(beta = safe_num(coef(fit)[2]), pred = predict(fit, newdata = test_df)))
    }
    if (method == "Theil-Sen") {
      fit <- mblm::mblm(log_y ~ log_x, data = train_df, repeated = TRUE)
      beta <- safe_num(coef(fit)[2])
      intercept <- safe_num(coef(fit)[1])
      return(list(beta = beta, pred = intercept + beta * test_df$log_x))
    }
    if (method == "Huber") {
      fit <- MASS::rlm(log_y ~ log_x, data = train_df, psi = MASS::psi.huber, maxit = 100)
      return(list(beta = safe_num(coef(fit)[2]), pred = predict(fit, newdata = test_df)))
    }
    if (method == "Quantile") {
      fit <- quantreg::rq(log_y ~ log_x, data = train_df, tau = 0.5)
      return(list(beta = safe_num(coef(fit)[2]), pred = predict(fit, newdata = test_df)))
    }
    list(beta = NA_real_, pred = rep(NA_real_, nrow(test_df)))
  }, error = function(e) list(beta = NA_real_, pred = rep(NA_real_, nrow(test_df))))
}

ensemble_slope <- function(df, k = 5) {
  n <- nrow(df)
  k <- min(k, n)
  folds <- sample(rep(seq_len(k), length.out = n))
  components <- c("OLS", "Theil-Sen", "Huber", "Quantile")
  cv_rmse <- rep(NA_real_, length(components))
  names(cv_rmse) <- components

  for (m in components) {
    pred_all <- rep(NA_real_, n)
    for (fold in seq_len(k)) {
      train_df <- df[folds != fold, , drop = FALSE]
      test_df  <- df[folds == fold, , drop = FALSE]
      out <- fit_predict_component(m, train_df, test_df)
      pred_all[folds == fold] <- out$pred
    }
    cv_rmse[m] <- sqrt(mean((df$log_y - pred_all)^2, na.rm = TRUE))
  }

  full_betas <- sapply(components, function(m) fit_predict_component(m, df, df)$beta)
  ok <- is.finite(cv_rmse) & cv_rmse > 0 & is.finite(full_betas)
  if (sum(ok) < 2) return(NA_real_)

  w <- 1 / (cv_rmse[ok]^2)
  w <- w / sum(w)
  sum(w * full_betas[ok])
}

# ----------------------------------------------------------------------------
# 4. SINGLE-HABITAT DATA-GENERATING PROCESS
# ----------------------------------------------------------------------------

simulate_single_habitat <- function(beta_true, n_obs, sigma_x_me, sigma_y_me, sim_id) {
  set.seed(1000000 + sim_id + round(1000 * beta_true) + 10 * n_obs +
             round(100 * sigma_x_me) + round(100 * sigma_y_me))

  log_x_true <- runif(n_obs, log(x_min), log(x_max))
  log_y_true <- log(alpha_true) + beta_true * log_x_true + rnorm(n_obs, 0, sigma_proc)

  log_x <- log_x_true + rnorm(n_obs, 0, sigma_x_me)
  log_y <- log_y_true + rnorm(n_obs, 0, sigma_y_me)

  tibble(
    log_x_true = log_x_true,
    log_y_true = log_y_true,
    log_x = log_x,
    log_y = log_y,
    x_obs = exp(log_x),
    y_obs = exp(log_y),
    habitat = factor("single")
  )
}

# ----------------------------------------------------------------------------
# 5. FIT EXACT MANUSCRIPT SINGLE-HABITAT METHODS
# ----------------------------------------------------------------------------

fit_single_methods <- function(df, beta_true, n_obs, sigma_x_me, sigma_y_me, sim_id) {
  x <- df$log_x
  y <- df$log_y

  out <- as.list(rep(NA_real_, length(METHOD_COLS)))
  names(out) <- names(METHOD_COLS)

  # OLS
  ols_fit <- tryCatch(lm(log_y ~ log_x, data = df), error = function(e) NULL)
  if (!is.null(ols_fit)) {
    out$beta_ols <- safe_num(coef(ols_fit)[2])
    # Clustered SE is a standard-error correction, not a different slope estimator.
    # In the single-habitat benchmark it shares the same beta estimate as OLS.
    out$beta_clustered_se <- out$beta_ols
    # With one habitat, random-intercept/slope models collapse to the global slope target.
    # The manuscript includes these rows for class comparison, so they are reported as OLS-equivalent.
    out$beta_lmm_random_intercept <- out$beta_ols
    out$beta_lmm_random_slope <- out$beta_ols
  }

  # ODR / weighted orthogonal regression using known axis-error SDs
  out$beta_odr <- odr_known_error_slope(x, y, sigma_x_me, sigma_y_me)

  # SMA Major Axis, RMA via lmodel2
  tryCatch({
    lm2 <- suppressMessages(lmodel2::lmodel2(log_y ~ log_x, data = df,
                                             range.y = "interval", range.x = "interval", nperm = 0))
    reg <- lm2$regression.results
    if (!is.null(reg)) {
      ma_idx <- which(tolower(reg$Method) %in% c("major axis", "ma"))[1]
      rma_idx <- which(grepl("reduced", tolower(reg$Method)) | tolower(reg$Method) == "rma")[1]
      if (!is.na(ma_idx)) out$beta_sma_major_axis <- safe_num(reg$Slope[ma_idx])
      if (!is.na(rma_idx)) out$beta_rma <- safe_num(reg$Slope[rma_idx])
    }
  }, error = function(e) {})

  # SMA via smatr
  tryCatch({
    sm <- smatr::sma(log_y ~ log_x, data = df)
    if (!is.null(sm$coef)) {
      cf <- sm$coef[[1]]
      out$beta_sma <- safe_num(cf["slope"])
    }
  }, error = function(e) {})
  if (is.na(out$beta_sma)) out$beta_sma <- sign(stats::cor(x, y)) * sd(y) / sd(x)

  # Deming
  out$beta_deming <- deming_slope(x, y, sigma_x_me, sigma_y_me)

  # Passing-Bablok
  tryCatch({
    pb <- mcr::mcreg(x, y, method.reg = "PaBa")
    out$beta_passing_bablok <- safe_num(pb@para[2, 1])
  }, error = function(e) {})

  # Theil-Sen, as described in the manuscript implementation
  tryCatch({
    ts <- mblm::mblm(log_y ~ log_x, data = df, repeated = TRUE)
    out$beta_theilsen <- safe_num(coef(ts)[2])
  }, error = function(e) {})

  # Median and trimmed pairwise slopes
  out$beta_median <- median_slope(x, y)
  out$beta_trimmed_mean <- trimmed_slope(x, y, trim = 0.10)

  # Huber
  tryCatch({
    hb <- MASS::rlm(log_y ~ log_x, data = df, psi = MASS::psi.huber, maxit = 100)
    out$beta_huber <- safe_num(coef(hb)[2])
  }, error = function(e) {})

  # Median quantile regression
  tryCatch({
    rq_fit <- quantreg::rq(log_y ~ log_x, data = df, tau = 0.5)
    out$beta_quantile <- safe_num(coef(rq_fit)[2])
  }, error = function(e) {})

  # Ridge / Lasso / Elastic Net / Adaptive Lasso
  tryCatch({
    x_mat <- as.matrix(df$log_x)
    y_vec <- df$log_y
    colnames(x_mat) <- "log_x"

    ridge <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 0, nfolds = min(5, nrow(df)), standardize = TRUE)
    lasso <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 1, nfolds = min(5, nrow(df)), standardize = TRUE)
    enet  <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 0.5, nfolds = min(5, nrow(df)), standardize = TRUE)

    out$beta_ridge <- safe_num(as.matrix(coef(ridge, s = "lambda.min"))["log_x", 1])
    out$beta_lasso <- safe_num(as.matrix(coef(lasso, s = "lambda.min"))["log_x", 1])
    out$beta_elastic_net <- safe_num(as.matrix(coef(enet, s = "lambda.min"))["log_x", 1])

    ols_b <- ifelse(is.finite(out$beta_ols) && abs(out$beta_ols) > 1e-6, abs(out$beta_ols), 1)
    penalty_factor <- 1 / ols_b
    adlasso <- glmnet::cv.glmnet(x_mat, y_vec, alpha = 1, nfolds = min(5, nrow(df)),
                                 standardize = TRUE, penalty.factor = penalty_factor)
    out$beta_adaptive_lasso <- safe_num(as.matrix(coef(adlasso, s = "lambda.min"))["log_x", 1])
  }, error = function(e) {})

  # SVM linear and radial basis function; beta is extracted from prediction-grid slope
  tryCatch({
    svm_lin <- e1071::svm(log_y ~ log_x, data = df, kernel = "linear", scale = TRUE)
    out$beta_svm_linear <- prediction_grid_slope(
      svm_lin,
      function(model, nd) predict(model, newdata = nd),
      min(x), max(x)
    )
  }, error = function(e) {})

  tryCatch({
    svm_rbf <- e1071::svm(log_y ~ log_x, data = df, kernel = "radial", scale = TRUE)
    out$beta_svm_rbf <- prediction_grid_slope(
      svm_rbf,
      function(model, nd) predict(model, newdata = nd),
      min(x), max(x)
    )
  }, error = function(e) {})

  # Ensemble: OLS + Theil-Sen + Huber + Quantile with CV inverse-RMSE weights
  out$beta_ensemble <- ensemble_slope(df, k = 5)

  as_tibble(c(
    list(
      sim_id = sim_id,
      beta_true = beta_true,
      n_obs = n_obs,
      sigma_x_me = sigma_x_me,
      sigma_y_me = sigma_y_me,
      error_ratio = sigma_y_me / sigma_x_me,
      total_error = sqrt(sigma_x_me^2 + sigma_y_me^2)
    ),
    out
  ))
}

# ----------------------------------------------------------------------------
# 6. SUMMARIES AND FIGURES
# ----------------------------------------------------------------------------

summarise_performance <- function(results_df) {
  results_df %>%
    pivot_longer(cols = all_of(names(METHOD_COLS)), names_to = "method_col", values_to = "estimate") %>%
    mutate(method = unname(METHOD_COLS[method_col])) %>%
    filter(!is.na(estimate), is.finite(estimate)) %>%
    group_by(n_obs, sigma_x_me, sigma_y_me, error_ratio, beta_true, method) %>%
    summarise(
      n_success = n(),
      mean_estimate = mean(estimate),
      bias = mean(estimate - beta_true),
      abs_bias = abs(mean(estimate - beta_true)),
      mae = mean(abs(estimate - beta_true)),
      rmse = sqrt(mean((estimate - beta_true)^2)),
      sd_estimate = sd(estimate),
      .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = METHOD_ORDER)) %>%
    arrange(beta_true, method)
}

summarise_overall <- function(perf_df) {
  perf_df %>%
    group_by(n_obs, sigma_x_me, sigma_y_me, error_ratio, method) %>%
    summarise(
      mean_abs_bias = mean(abs_bias, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      score = mean_abs_bias + mean_rmse,
      .groups = "drop"
    ) %>%
    arrange(score) %>%
    mutate(rank = row_number())
}

save_figures <- function(perf_df, overall_df, scenario_dir, n_obs, sigma_x_me, sigma_y_me) {
  p_rmse <- ggplot(perf_df, aes(x = factor(beta_true), y = method, fill = rmse)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.3f", rmse)), size = 2.6) +
    scale_y_discrete(limits = rev(METHOD_ORDER), drop = FALSE) +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), panel.grid = element_blank()) +
    labs(
      title = sprintf("Single-habitat RMSE: n=%d, sigma_x=%.2f, sigma_y=%.2f", n_obs, sigma_x_me, sigma_y_me),
      x = "True beta", y = "Method", fill = "RMSE"
    )

  p_bias <- ggplot(perf_df, aes(x = factor(beta_true), y = method, fill = bias)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.3f", bias)), size = 2.6) +
    scale_y_discrete(limits = rev(METHOD_ORDER), drop = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), panel.grid = element_blank()) +
    labs(
      title = sprintf("Single-habitat bias: n=%d, sigma_x=%.2f, sigma_y=%.2f", n_obs, sigma_x_me, sigma_y_me),
      x = "True beta", y = "Method", fill = "Bias"
    )

  p_rank <- ggplot(overall_df, aes(x = reorder(method, mean_rmse), y = mean_rmse)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      title = sprintf("Single-habitat overall ranking: n=%d, sigma_x=%.2f, sigma_y=%.2f", n_obs, sigma_x_me, sigma_y_me),
      x = "Method", y = "Mean RMSE across beta"
    )

  ggsave(file.path(scenario_dir, "single_habitat_rmse_heatmap.png"), p_rmse, width = 8.5, height = 8.5, dpi = 300)
  ggsave(file.path(scenario_dir, "single_habitat_bias_heatmap.png"), p_bias, width = 8.5, height = 8.5, dpi = 300)
  ggsave(file.path(scenario_dir, "single_habitat_method_ranking.png"), p_rank, width = 8, height = 7, dpi = 300)
}

# ----------------------------------------------------------------------------
# 7. RUN SCENARIOS AND SAVE AFTER EACH REPLICATE
# ----------------------------------------------------------------------------

all_combined <- list()
scenario_counter <- 0
n_scenarios <- length(n_values) * nrow(error_configs)

cat("============================================================\n")
cat("SINGLE-HABITAT MANUSCRIPT BENCHMARK\n")
cat("Exact single-habitat method rows, no extra method classes\n")
cat("============================================================\n")
cat(sprintf("Replicates per beta/scenario: %d\n", n_sim))
cat(sprintf("Output base folder: %s\n\n", base_out_dir))

for (n_obs in n_values) {
  for (ec_i in seq_len(nrow(error_configs))) {

    sigma_x_me <- error_configs$sigma_x_me[ec_i]
    sigma_y_me <- error_configs$sigma_y_me[ec_i]

    scenario_counter <- scenario_counter + 1
    this_folder <- folder_name(n_obs, sigma_x_me, sigma_y_me)
    scenario_dir <- file.path(base_out_dir, this_folder)
    if (!dir.exists(scenario_dir)) dir.create(scenario_dir, recursive = TRUE)

    cat("------------------------------------------------------------\n")
    cat(sprintf("Scenario %d/%d: %s\n", scenario_counter, n_scenarios, this_folder))
    cat("------------------------------------------------------------\n")

    running_file <- file.path(scenario_dir, "all_results_running.csv")
    if (file.exists(running_file)) file.remove(running_file)

    scenario_results <- list()
    row_i <- 0

    for (beta_true in beta_values) {
      cat(sprintf("  beta = %.1f\n", beta_true))

      for (sim_id in seq_len(n_sim)) {
        df <- simulate_single_habitat(
          beta_true = beta_true,
          n_obs = n_obs,
          sigma_x_me = sigma_x_me,
          sigma_y_me = sigma_y_me,
          sim_id = sim_id
        )

        res <- fit_single_methods(
          df = df,
          beta_true = beta_true,
          n_obs = n_obs,
          sigma_x_me = sigma_x_me,
          sigma_y_me = sigma_y_me,
          sim_id = sim_id
        )

        row_i <- row_i + 1
        scenario_results[[row_i]] <- res

        # Save after each replicate to avoid losing progress.
        readr::write_csv(bind_rows(scenario_results), running_file)

        if (sim_id %% 25 == 0 || sim_id == n_sim) {
          cat(sprintf("    completed %d/%d\n", sim_id, n_sim))
        }
      }
    }

    scenario_df <- bind_rows(scenario_results)
    scenario_perf <- summarise_performance(scenario_df)
    scenario_overall <- summarise_overall(scenario_perf)

    readr::write_csv(scenario_df, file.path(scenario_dir, "all_results.csv"))
    readr::write_csv(scenario_perf, file.path(scenario_dir, "performance_by_beta.csv"))
    readr::write_csv(scenario_overall, file.path(scenario_dir, "overall_ranking.csv"))
    saveRDS(
      list(
        raw = scenario_df,
        performance_by_beta = scenario_perf,
        overall_ranking = scenario_overall,
        params = list(n_sim = n_sim, n_obs = n_obs, beta_values = beta_values,
                      sigma_x_me = sigma_x_me, sigma_y_me = sigma_y_me,
                      alpha_true = alpha_true, sigma_proc = sigma_proc)
      ),
      file.path(scenario_dir, "complete_results.rds")
    )

    save_figures(scenario_perf, scenario_overall, scenario_dir, n_obs, sigma_x_me, sigma_y_me)

    all_combined[[length(all_combined) + 1]] <- scenario_df
  }
}

# ----------------------------------------------------------------------------
# 8. COMBINED OUTPUTS ACROSS ALL SIX FOLDERS
# ----------------------------------------------------------------------------

combined_results <- bind_rows(all_combined)
combined_performance <- summarise_performance(combined_results)
combined_overall <- combined_performance %>%
  group_by(method) %>%
  summarise(
    AbsBias = mean(abs_bias, na.rm = TRUE),
    RMSE = mean(rmse, na.rm = TRUE),
    Score = AbsBias + RMSE,
    .groups = "drop"
  ) %>%
  arrange(Score) %>%
  mutate(rank = row_number())

readr::write_csv(combined_results, file.path(base_out_dir, "combined_all_results.csv"))
readr::write_csv(combined_performance, file.path(base_out_dir, "combined_performance_by_beta_n_error.csv"))
readr::write_csv(combined_overall, file.path(base_out_dir, "combined_overall_ranking.csv"))
saveRDS(
  list(
    raw = combined_results,
    performance = combined_performance,
    overall = combined_overall,
    method_order = METHOD_ORDER,
    params = list(n_sim = n_sim, n_values = n_values, beta_values = beta_values,
                  error_configs = error_configs, alpha_true = alpha_true, sigma_proc = sigma_proc)
  ),
  file.path(base_out_dir, "single_habitat_exact_methods_complete.rds")
)

# Combined heatmap matching supplementary-style scenario display
combined_plot_df <- combined_performance %>%
  mutate(
    scenario = sprintf("beta=%s_n=%d_eta=%s", beta_true, n_obs, error_ratio),
    method = factor(method, levels = rev(METHOD_ORDER))
  )

p_comb_rmse <- ggplot(combined_plot_df, aes(x = scenario, y = method, fill = rmse)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7), panel.grid = element_blank()) +
  labs(title = "Single-habitat scenario-specific RMSE heatmap", x = "Scenario", y = "Method", fill = "RMSE")

ggsave(file.path(base_out_dir, "combined_single_habitat_rmse_heatmap.png"),
       p_comb_rmse, width = 13, height = 8, dpi = 300)

cat("\n============================================================\n")
cat("COMPLETE: SINGLE-HABITAT MANUSCRIPT BENCHMARK\n")
cat("============================================================\n")
cat(sprintf("All folders saved under: %s\n", base_out_dir))
cat("Folders created:\n")
print(list.dirs(base_out_dir, recursive = FALSE, full.names = FALSE))
cat("\nCombined ranking:\n")
print(combined_overall)
cat("============================================================\n")
