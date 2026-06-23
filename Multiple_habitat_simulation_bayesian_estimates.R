
# ============================================================
# Targeted multi-habitat extension:
# Top-performing manuscript models + Bayesian hierarchical EIV
#
# Purpose:
#   This script does NOT replace the main 20-method benchmark.
#   It provides a focused follow-up analysis for the Discussion/
#   Future Directions section:
#
#     "Can a fully Bayesian hierarchical measurement-error model
#      improve estimation when habitat hierarchy and EIV coexist?"
#
# Models compared:
#   1. OLS
#   2. Theil-Sen
#   3. Siegel
#   4. Huber
#   5. Quantile regression
#   6. Ensemble
#   7. MixedLM without EIV
#   8. Deming/global EIV
#   9. Bayesian hierarchical mixed-effects EIV model via CmdStanR
#
# No PyStan. No CmdStanPy.
# Bayesian model is fitted using R + CmdStanR.
# ============================================================


# ============================================================
# 0. PACKAGE SETUP
# ============================================================

# Install once if needed:
# install.packages(
#   c("cmdstanr", "lme4", "MASS", "dplyr", "readr", "posterior", "quantreg"),
#   repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
# )

suppressPackageStartupMessages({
  library(cmdstanr)
  library(lme4)
  library(MASS)
  library(dplyr)
  library(readr)
  library(posterior)
  library(quantreg)
})

options(cmdstanr_no_ver_check = TRUE)


# ============================================================
# 1. CMDSTAN SETUP
# ============================================================

setup_cmdstan_clean <- function() {

  # Prefer clean CmdStanR installations over conda CmdStan.
  candidate_dirs <- character(0)

  custom_parent <- "C:/Users/hanra/cmdstanr_cmdstan"
  if (dir.exists(custom_parent)) {
    custom_dirs <- list.dirs(custom_parent, recursive = FALSE, full.names = TRUE)
    custom_dirs <- custom_dirs[grepl("cmdstan-", basename(custom_dirs))]
    candidate_dirs <- c(candidate_dirs, custom_dirs)
  }

  default_parent <- file.path(Sys.getenv("USERPROFILE"), ".cmdstan")
  if (dir.exists(default_parent)) {
    default_dirs <- list.dirs(default_parent, recursive = FALSE, full.names = TRUE)
    default_dirs <- default_dirs[grepl("cmdstan-", basename(default_dirs))]
    candidate_dirs <- c(candidate_dirs, default_dirs)
  }

  candidate_dirs <- unique(candidate_dirs)
  candidate_dirs <- candidate_dirs[dir.exists(candidate_dirs)]

  # Avoid the conda CmdStan path that previously caused the print-STANCFLAGS error.
  candidate_dirs <- candidate_dirs[
    !grepl("anaconda3/envs/stan_env", candidate_dirs, ignore.case = TRUE)
  ]
  candidate_dirs <- candidate_dirs[
    !grepl("anaconda3\\\\envs\\\\stan_env", candidate_dirs, ignore.case = TRUE)
  ]

  if (length(candidate_dirs) == 0) {
    stop(
      "No clean CmdStan installation found.\n",
      "Run this first in R:\n",
      "dir.create('C:/Users/hanra/cmdstanr_cmdstan', recursive = TRUE, showWarnings = FALSE)\n",
      "cmdstanr::install_cmdstan(dir = 'C:/Users/hanra/cmdstanr_cmdstan', overwrite = TRUE)\n"
    )
  }

  candidate_dirs <- sort(candidate_dirs)
  selected <- candidate_dirs[length(candidate_dirs)]

  cmdstanr::set_cmdstan_path(selected)

  cat("\nUsing CmdStan path:\n")
  print(cmdstanr::cmdstan_path())

  cat("\nCmdStan version:\n")
  print(cmdstanr::cmdstan_version())

  cat("\nChecking C++ toolchain:\n")
  cmdstanr::check_cmdstan_toolchain()

  # Safety test.
  cat("\nTesting CmdStan with Bernoulli example...\n")

  bernoulli_file <- file.path(
    cmdstanr::cmdstan_path(),
    "examples",
    "bernoulli",
    "bernoulli.stan"
  )

  if (!file.exists(bernoulli_file)) {
    stop("Bernoulli example not found at: ", bernoulli_file)
  }

  invisible(
    tryCatch(
      cmdstanr::cmdstan_model(bernoulli_file, force_recompile = FALSE),
      error = function(e) {
        stop(
          "CmdStan Bernoulli test failed. Do not continue.\n",
          "Error was:\n",
          conditionMessage(e)
        )
      }
    )
  )

  cat("Bernoulli test compiled successfully.\n\n")
}

setup_cmdstan_clean()


# ============================================================
# 2. GLOBAL SETTINGS
# ============================================================

OUT_DIR <- "targeted_multihabitat_bayesian_eiv_extension"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

set.seed(123)

# Start small. Increase after confirming runtime.
# Suggested:
#   quick test: 3
#   preliminary section: 10 or 20
#   manuscript-quality extension: 50 or 100
N_REPLICATES <- 5

TRUE_BETAS <- c(0.6, 1.0, 1.4)

ERROR_SETTINGS <- data.frame(
  sigma_x = c(0.20, 0.30, 0.40),
  sigma_y = c(0.40, 0.30, 0.20),
  label   = c("eta_2.0", "eta_1.0", "eta_0.5")
)

# Manuscript-style multi-habitat scenarios:
# H = 6 or 12; total n = 36 or 72.
H_VALUES <- c(6, 12)
N_TOTAL_VALUES <- c(36, 72)

ALPHA0_TRUE <- 0.5
SIGMA_EPS_TRUE <- 0.10

SIGMA_ALPHA_TRUE <- 0.40
SIGMA_BETA_TRUE <- 0.40
RHO_AB_TRUE <- -0.30

SIGMA_XSTAR_TRUE <- 0.40

HABITAT_MEANS_6  <- c(1.0, 2.5, 4.0, 5.5, 7.0, 8.5)
HABITAT_MEANS_12 <- c(1.0, 1.8, 2.6, 3.4, 4.2, 5.0,
                      5.8, 6.6, 7.4, 8.2, 9.0, 9.8)

# Bayesian sampling settings.
# Increase later if needed.
STAN_CHAINS <- 4
STAN_WARMUP <- 500
STAN_SAMPLING <- 500

# Ensemble settings.
ENSEMBLE_KFOLD <- 5


# ============================================================
# 3. STAN MODEL CODE: BAYESIAN HIERARCHICAL EIV
# ============================================================

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> H;
  array[N] int<lower=1, upper=H> habitat;

  vector[N] x_obs;
  vector[N] y_obs;

  vector[H] mu_x_h;

  real<lower=0> sigma_x_meas;
  real<lower=0> sigma_y_meas;

  real<lower=0> beta_bound;
  real<lower=0> alpha0_bound;
  real<lower=0> re_sd_bound;
  real<lower=0> sigma_eps_bound;
  real<lower=0> sigma_xstar_bound;
}

parameters {
  real<lower=-alpha0_bound, upper=alpha0_bound> alpha0;
  real<lower=-beta_bound, upper=beta_bound> beta;

  vector[N] x_star;

  vector[H] alpha_raw;
  vector[H] beta_raw;

  real<lower=0, upper=re_sd_bound> sigma_alpha;
  real<lower=0, upper=re_sd_bound> sigma_beta;

  real<lower=-0.95, upper=0.95> rho_ab;

  real<lower=0, upper=sigma_eps_bound> sigma_eps;
  real<lower=0, upper=sigma_xstar_bound> sigma_xstar;
}

transformed parameters {
  vector[H] alpha_h;
  vector[H] beta_h;

  alpha_h = sigma_alpha * alpha_raw;

  beta_h = sigma_beta * (
    rho_ab * alpha_raw +
    sqrt(1 - square(rho_ab)) * beta_raw
  );
}

model {
  // Broad uniform priors are implemented through parameter bounds
  // for alpha0, beta, variance components, rho_ab, and sigma_xstar.

  // These are weak regularizing priors for non-centered standardized
  // random-effect coordinates only.
  alpha_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);

  for (n in 1:N) {
    x_star[n] ~ normal(mu_x_h[habitat[n]], sigma_xstar);
  }

  x_obs ~ normal(x_star, sigma_x_meas);

  for (n in 1:N) {
    y_obs[n] ~ normal(
      alpha0 +
      alpha_h[habitat[n]] +
      (beta + beta_h[habitat[n]]) * x_star[n],
      sqrt(square(sigma_y_meas) + square(sigma_eps))
    );
  }
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  for (n in 1:N) {
    real mu_y;
    real sigma_y_total;

    mu_y =
      alpha0 +
      alpha_h[habitat[n]] +
      (beta + beta_h[habitat[n]]) * x_star[n];

    sigma_y_total = sqrt(square(sigma_y_meas) + square(sigma_eps));

    y_rep[n] = normal_rng(mu_y, sigma_y_total);
    log_lik[n] = normal_lpdf(y_obs[n] | mu_y, sigma_y_total);
  }
}
"

stan_file <- file.path(OUT_DIR, "targeted_multihabitat_bayes_eiv.stan")
writeLines(stan_code, stan_file)


# ============================================================
# 4. MULTI-HABITAT SIMULATION
# ============================================================

get_habitat_means <- function(H) {
  if (H == 6) return(HABITAT_MEANS_6)
  if (H == 12) return(HABITAT_MEANS_12)
  stop("Unsupported H: ", H)
}

simulate_multihabitat_data <- function(
    beta_true,
    sigma_x,
    sigma_y,
    H,
    n_per,
    alpha0 = ALPHA0_TRUE,
    sigma_eps = SIGMA_EPS_TRUE,
    sigma_alpha = SIGMA_ALPHA_TRUE,
    sigma_beta = SIGMA_BETA_TRUE,
    rho = RHO_AB_TRUE,
    sigma_xstar = SIGMA_XSTAR_TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  habitat_means <- get_habitat_means(H)

  cov_re <- matrix(
    c(
      sigma_alpha^2,
      rho * sigma_alpha * sigma_beta,
      rho * sigma_alpha * sigma_beta,
      sigma_beta^2
    ),
    nrow = 2,
    byrow = TRUE
  )

  random_effects <- MASS::mvrnorm(
    n = H,
    mu = c(0, 0),
    Sigma = cov_re
  )

  rows <- list()
  counter <- 1

  for (h in seq_len(H)) {

    alpha_h <- random_effects[h, 1]
    beta_h  <- random_effects[h, 2]
    mu_h    <- habitat_means[h]

    x_star <- rnorm(n_per, mean = mu_h, sd = sigma_xstar)
    eps    <- rnorm(n_per, mean = 0, sd = sigma_eps)

    y_star <- alpha0 + alpha_h + (beta_true + beta_h) * x_star + eps

    x_obs <- x_star + rnorm(n_per, mean = 0, sd = sigma_x)
    y_obs <- y_star + rnorm(n_per, mean = 0, sd = sigma_y)

    for (j in seq_len(n_per)) {
      rows[[counter]] <- data.frame(
        habitat = h,
        j = j,
        x_star = x_star[j],
        y_star = y_star[j],
        x_obs = x_obs[j],
        y_obs = y_obs[j],
        beta_true = beta_true,
        alpha_h_true = alpha_h,
        beta_h_true = beta_h
      )
      counter <- counter + 1
    }
  }

  dplyr::bind_rows(rows)
}


# ============================================================
# 5. HELPER FUNCTIONS FOR ROBUST MODELS
# ============================================================

safe_slope_intercept <- function(beta, x, y) {
  alpha <- median(y - beta * x, na.rm = TRUE)
  c(alpha = as.numeric(alpha), beta = as.numeric(beta))
}

estimate_ols_pair <- function(df) {
  fit <- lm(y_obs ~ x_obs, data = df)
  c(alpha = as.numeric(coef(fit)[1]), beta = as.numeric(coef(fit)[2]))
}

predict_pair <- function(coef_pair, newdata) {
  as.numeric(coef_pair["alpha"] + coef_pair["beta"] * newdata$x_obs)
}

pairwise_slopes <- function(x, y) {
  n <- length(x)
  if (n < 2) return(numeric(0))

  pairs <- combn(seq_len(n), 2)
  dx <- x[pairs[2, ]] - x[pairs[1, ]]
  dy <- y[pairs[2, ]] - y[pairs[1, ]]

  ok <- is.finite(dx) & is.finite(dy) & abs(dx) > 1e-12
  dy[ok] / dx[ok]
}

estimate_theilsen_pair <- function(df) {
  x <- df$x_obs
  y <- df$y_obs
  slopes <- pairwise_slopes(x, y)

  if (length(slopes) == 0) return(c(alpha = NA_real_, beta = NA_real_))

  beta <- median(slopes, na.rm = TRUE)
  safe_slope_intercept(beta, x, y)
}

estimate_siegel_pair <- function(df) {
  x <- df$x_obs
  y <- df$y_obs
  n <- length(x)

  med_slopes <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    dx <- x[-i] - x[i]
    dy <- y[-i] - y[i]
    ok <- is.finite(dx) & is.finite(dy) & abs(dx) > 1e-12

    if (sum(ok) > 0) {
      med_slopes[i] <- median(dy[ok] / dx[ok], na.rm = TRUE)
    }
  }

  beta <- median(med_slopes, na.rm = TRUE)

  if (!is.finite(beta)) return(c(alpha = NA_real_, beta = NA_real_))

  safe_slope_intercept(beta, x, y)
}

estimate_huber_pair <- function(df) {
  fit <- tryCatch(
    MASS::rlm(y_obs ~ x_obs, data = df, psi = MASS::psi.huber, maxit = 100),
    error = function(e) NULL
  )

  if (is.null(fit)) return(c(alpha = NA_real_, beta = NA_real_))

  c(alpha = as.numeric(coef(fit)[1]), beta = as.numeric(coef(fit)[2]))
}

estimate_quantile_pair <- function(df) {
  fit <- tryCatch(
    quantreg::rq(y_obs ~ x_obs, data = df, tau = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) return(c(alpha = NA_real_, beta = NA_real_))

  c(alpha = as.numeric(coef(fit)[1]), beta = as.numeric(coef(fit)[2]))
}


# ============================================================
# 6. CORE ESTIMATORS
# ============================================================

estimate_ols <- function(df) {
  estimate_ols_pair(df)["beta"]
}

estimate_theilsen <- function(df) {
  estimate_theilsen_pair(df)["beta"]
}

estimate_siegel <- function(df) {
  estimate_siegel_pair(df)["beta"]
}

estimate_huber <- function(df) {
  estimate_huber_pair(df)["beta"]
}

estimate_quantile <- function(df) {
  estimate_quantile_pair(df)["beta"]
}

estimate_deming_global <- function(df, sigma_x, sigma_y) {

  x <- df$x_obs
  y <- df$y_obs

  sx2 <- var(x)
  sy2 <- var(y)
  sxy <- cov(x, y)

  if (!is.finite(sxy) || abs(sxy) < 1e-12) return(NA_real_)

  lambda <- sigma_y^2 / sigma_x^2

  beta_hat <- (
    sy2 - lambda * sx2 +
      sqrt((sy2 - lambda * sx2)^2 + 4 * lambda * sxy^2)
  ) / (2 * sxy)

  as.numeric(beta_hat)
}

estimate_mixedlm <- function(df) {

  df$habitat <- factor(df$habitat)

  fit <- tryCatch(
    {
      lme4::lmer(
        y_obs ~ x_obs + (1 + x_obs | habitat),
        data = df,
        REML = FALSE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 100000)
        )
      )
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    fit <- tryCatch(
      {
        lme4::lmer(
          y_obs ~ x_obs + (1 | habitat),
          data = df,
          REML = FALSE,
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 100000)
          )
        )
      },
      error = function(e) NULL
    )
  }

  if (is.null(fit)) return(NA_real_)

  as.numeric(lme4::fixef(fit)["x_obs"])
}


# ============================================================
# 7. ENSEMBLE ESTIMATOR
# ============================================================

fitters_for_ensemble <- list(
  OLS = estimate_ols_pair,
  Theil_Sen = estimate_theilsen_pair,
  Siegel = estimate_siegel_pair,
  Huber = estimate_huber_pair,
  Quantile = estimate_quantile_pair
)

estimate_ensemble <- function(df, kfold = ENSEMBLE_KFOLD, seed = 123) {

  set.seed(seed)

  n <- nrow(df)
  kfold <- min(kfold, n)

  fold_id <- sample(rep(seq_len(kfold), length.out = n))

  rmse_by_method <- rep(NA_real_, length(fitters_for_ensemble))
  names(rmse_by_method) <- names(fitters_for_ensemble)

  for (m in names(fitters_for_ensemble)) {

    pred_all <- rep(NA_real_, n)

    for (k in seq_len(kfold)) {

      train <- df[fold_id != k, , drop = FALSE]
      test  <- df[fold_id == k, , drop = FALSE]

      coef_pair <- tryCatch(
        fitters_for_ensemble[[m]](train),
        error = function(e) c(alpha = NA_real_, beta = NA_real_)
      )

      if (all(is.finite(coef_pair))) {
        pred_all[fold_id == k] <- predict_pair(coef_pair, test)
      }
    }

    ok <- is.finite(pred_all) & is.finite(df$y_obs)

    if (sum(ok) > 3) {
      rmse_by_method[m] <- sqrt(mean((df$y_obs[ok] - pred_all[ok])^2))
    }
  }

  full_slopes <- sapply(
    fitters_for_ensemble,
    function(f) {
      out <- tryCatch(f(df), error = function(e) c(alpha = NA_real_, beta = NA_real_))
      as.numeric(out["beta"])
    }
  )

  ok <- is.finite(rmse_by_method) &
    is.finite(full_slopes) &
    rmse_by_method > 0

  if (sum(ok) == 0) return(NA_real_)

  weights <- 1 / (rmse_by_method[ok]^2)
  weights <- weights / sum(weights)

  beta_ens <- sum(weights * full_slopes[ok])

  as.numeric(beta_ens)
}


# ============================================================
# 8. BAYESIAN CMDSTANR EIV ESTIMATOR
# ============================================================

estimate_cmdstanr_mixed_eiv <- function(
    df,
    sigma_x,
    sigma_y,
    model,
    seed = 123,
    chains = STAN_CHAINS,
    iter_warmup = STAN_WARMUP,
    iter_sampling = STAN_SAMPLING
) {

  df <- df %>% arrange(habitat, j)

  habitat_index <- as.integer(factor(df$habitat))
  H <- length(unique(habitat_index))

  mu_x_h <- df %>%
    mutate(habitat_index = habitat_index) %>%
    group_by(habitat_index) %>%
    summarise(mu_x = mean(x_obs), .groups = "drop") %>%
    arrange(habitat_index) %>%
    pull(mu_x)

  stan_data <- list(
    N = nrow(df),
    H = H,
    habitat = habitat_index,
    x_obs = df$x_obs,
    y_obs = df$y_obs,
    mu_x_h = mu_x_h,
    sigma_x_meas = sigma_x,
    sigma_y_meas = sigma_y,
    beta_bound = 5.0,
    alpha0_bound = 20.0,
    re_sd_bound = 3.0,
    sigma_eps_bound = 3.0,
    sigma_xstar_bound = 3.0
  )

  ols_beta <- estimate_ols(df)
  alpha0_init <- mean(df$y_obs) - ols_beta * mean(df$x_obs)

  init_fun <- function() {
    list(
      alpha0 = max(min(alpha0_init, 19.5), -19.5),
      beta = max(min(ols_beta, 4.5), -4.5),
      x_star = df$x_obs,
      alpha_raw = rnorm(H, 0, 0.1),
      beta_raw = rnorm(H, 0, 0.1),
      sigma_alpha = 0.4,
      sigma_beta = 0.4,
      rho_ab = -0.1,
      sigma_eps = 0.1,
      sigma_xstar = 0.4
    )
  }

  fit <- tryCatch(
    {
      model$sample(
        data = stan_data,
        seed = seed,
        chains = chains,
        parallel_chains = min(chains, 2),
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        init = init_fun,
        refresh = 0,
        adapt_delta = 0.95,
        max_treedepth = 12
      )
    },
    error = function(e) {
      message("Stan failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(data.frame(
      beta_hat = NA_real_,
      rhat = NA_real_,
      ess_bulk = NA_real_,
      n_divergent = NA_real_
    ))
  }

  beta_draws <- posterior::as_draws_df(fit$draws("beta"))$beta
  beta_hat <- mean(beta_draws)

  summ <- fit$summary("beta")

  sampler_diag <- tryCatch(
    as.data.frame(fit$sampler_diagnostics(format = "df")),
    error = function(e) NULL
  )

  n_divergent <- NA_real_
  if (!is.null(sampler_diag) && "divergent__" %in% names(sampler_diag)) {
    n_divergent <- sum(sampler_diag$divergent__, na.rm = TRUE)
  }

  data.frame(
    beta_hat = beta_hat,
    rhat = summ$rhat,
    ess_bulk = summ$ess_bulk,
    n_divergent = n_divergent
  )
}


# ============================================================
# 9. SUMMARY FUNCTIONS
# ============================================================

summarize_accuracy <- function(results_df) {

  results_df %>%
    group_by(beta_true, H, n_total, n_per, error_setting, method) %>%
    summarise(
      mean_beta_hat = mean(beta_hat, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      abs_bias = abs(mean(error, na.rm = TRUE)),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs(error), na.rm = TRUE),
      sd_beta_hat = sd(beta_hat, na.rm = TRUE),
      n_success = sum(is.finite(beta_hat)),
      n_total_reps = n(),
      failure_rate = 1 - n_success / n_total_reps,
      .groups = "drop"
    ) %>%
    group_by(beta_true, H, n_total, n_per, error_setting) %>%
    mutate(rank_within_scenario = rank(rmse, ties.method = "min")) %>%
    ungroup() %>%
    arrange(beta_true, H, n_total, error_setting, rmse)
}

summarize_overall <- function(results_df) {

  results_df %>%
    group_by(method) %>%
    summarise(
      mean_beta_hat = mean(beta_hat, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      abs_bias = abs(mean(error, na.rm = TRUE)),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs(error), na.rm = TRUE),
      sd_beta_hat = sd(beta_hat, na.rm = TRUE),
      n_success = sum(is.finite(beta_hat)),
      n_total_reps = n(),
      failure_rate = 1 - n_success / n_total_reps,
      .groups = "drop"
    ) %>%
    arrange(rmse)
}

summarize_by_beta <- function(results_df) {

  results_df %>%
    group_by(beta_true, method) %>%
    summarise(
      mean_beta_hat = mean(beta_hat, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      abs_bias = abs(mean(error, na.rm = TRUE)),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs(error), na.rm = TRUE),
      n_success = sum(is.finite(beta_hat)),
      n_total_reps = n(),
      failure_rate = 1 - n_success / n_total_reps,
      .groups = "drop"
    ) %>%
    arrange(beta_true, rmse)
}

summarize_by_design <- function(results_df) {

  results_df %>%
    group_by(H, n_total, n_per, method) %>%
    summarise(
      bias = mean(error, na.rm = TRUE),
      abs_bias = abs(mean(error, na.rm = TRUE)),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mae = mean(abs(error), na.rm = TRUE),
      n_success = sum(is.finite(beta_hat)),
      n_total_reps = n(),
      failure_rate = 1 - n_success / n_total_reps,
      .groups = "drop"
    ) %>%
    arrange(H, n_total, rmse)
}


# ============================================================
# 10. COMPILE BAYESIAN MODEL
# ============================================================

cat("\nCompiling Bayesian hierarchical EIV Stan model...\n")

bayes_model <- tryCatch(
  {
    cmdstanr::cmdstan_model(
      stan_file,
      force_recompile = TRUE
    )
  },
  error = function(e) {
    stop(
      "Stan compilation failed. Do not continue.\n",
      "Error was:\n",
      conditionMessage(e)
    )
  }
)

cat("Bayesian EIV Stan model compiled successfully.\n\n")


# ============================================================
# 11. SCENARIO GRID
# ============================================================

scenario_grid <- expand.grid(
  beta_true = TRUE_BETAS,
  H = H_VALUES,
  n_total = N_TOTAL_VALUES,
  error_index = seq_len(nrow(ERROR_SETTINGS)),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    n_per = n_total / H,
    sigma_x = ERROR_SETTINGS$sigma_x[error_index],
    sigma_y = ERROR_SETTINGS$sigma_y[error_index],
    error_setting = ERROR_SETTINGS$label[error_index]
  ) %>%
  filter(abs(n_per - round(n_per)) < 1e-12) %>%
  mutate(n_per = as.integer(n_per)) %>%
  arrange(beta_true, H, n_total, error_setting)

cat("Scenario grid:\n")
print(scenario_grid)

cat("\nTotal scenarios:", nrow(scenario_grid), "\n")
cat("Replicates per scenario:", N_REPLICATES, "\n\n")


# ============================================================
# 12. MAIN LOOP
# ============================================================

all_rows <- list()
diag_rows <- list()
row_counter <- 1
diag_counter <- 1

cat("Running targeted multi-habitat extension...\n\n")

for (s in seq_len(nrow(scenario_grid))) {

  sc <- scenario_grid[s, ]

  beta_true <- sc$beta_true
  H <- sc$H
  n_total <- sc$n_total
  n_per <- sc$n_per
  sigma_x <- sc$sigma_x
  sigma_y <- sc$sigma_y
  label <- sc$error_setting

  for (r in seq_len(N_REPLICATES)) {

    seed <- 1000000 +
      s * 1000 +
      r +
      as.integer(beta_true * 100) +
      as.integer(H * 10) +
      as.integer(n_total)

    df <- simulate_multihabitat_data(
      beta_true = beta_true,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      H = H,
      n_per = n_per,
      seed = seed
    )

    beta_ols <- estimate_ols(df)
    beta_ts <- estimate_theilsen(df)
    beta_siegel <- estimate_siegel(df)
    beta_huber <- estimate_huber(df)
    beta_quantile <- estimate_quantile(df)
    beta_ens <- estimate_ensemble(df, kfold = ENSEMBLE_KFOLD, seed = seed + 33)
    beta_mixed <- estimate_mixedlm(df)
    beta_dem <- estimate_deming_global(df, sigma_x = sigma_x, sigma_y = sigma_y)

    stan_res <- estimate_cmdstanr_mixed_eiv(
      df = df,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      model = bayes_model,
      seed = seed,
      chains = STAN_CHAINS,
      iter_warmup = STAN_WARMUP,
      iter_sampling = STAN_SAMPLING
    )

    estimates <- data.frame(
      method = c(
        "OLS",
        "Theil_Sen",
        "Siegel",
        "Huber",
        "Quantile",
        "Ensemble",
        "MixedLM_no_EIV",
        "Deming_global_EIV",
        "Bayes_Mixed_EIV_CmdStanR"
      ),
      beta_hat = c(
        beta_ols,
        beta_ts,
        beta_siegel,
        beta_huber,
        beta_quantile,
        beta_ens,
        beta_mixed,
        beta_dem,
        stan_res$beta_hat
      )
    )

    for (k in seq_len(nrow(estimates))) {

      beta_hat <- estimates$beta_hat[k]

      all_rows[[row_counter]] <- data.frame(
        scenario_id = s,
        replicate = r,
        beta_true = beta_true,
        H = H,
        n_total = n_total,
        n_per = n_per,
        sigma_x = sigma_x,
        sigma_y = sigma_y,
        error_setting = label,
        method = estimates$method[k],
        beta_hat = beta_hat,
        error = ifelse(is.finite(beta_hat), beta_hat - beta_true, NA_real_)
      )

      row_counter <- row_counter + 1
    }

    diag_rows[[diag_counter]] <- data.frame(
      scenario_id = s,
      replicate = r,
      beta_true = beta_true,
      H = H,
      n_total = n_total,
      n_per = n_per,
      error_setting = label,
      stan_beta_rhat = stan_res$rhat,
      stan_beta_ess_bulk = stan_res$ess_bulk,
      stan_n_divergent = stan_res$n_divergent
    )

    diag_counter <- diag_counter + 1

    cat(
      "Scenario ", s, "/", nrow(scenario_grid),
      " | beta=", beta_true,
      " | H=", H,
      " | n=", n_total,
      " | nper=", n_per,
      " | ", label,
      " | rep=", r, "/", N_REPLICATES,
      " | Bayes beta=", round(stan_res$beta_hat, 3),
      " | Rhat=", round(stan_res$rhat, 3),
      " | ESS=", round(stan_res$ess_bulk, 1),
      " | div=", stan_res$n_divergent,
      "\n",
      sep = ""
    )
  }
}


# ============================================================
# 13. SAVE RESULTS
# ============================================================

results_df <- dplyr::bind_rows(all_rows)
diag_df <- dplyr::bind_rows(diag_rows)

scenario_summary <- summarize_accuracy(results_df)
overall_summary <- summarize_overall(results_df)
beta_summary <- summarize_by_beta(results_df)
design_summary <- summarize_by_design(results_df)

results_path <- file.path(OUT_DIR, "targeted_multihabitat_detailed_results.csv")
scenario_path <- file.path(OUT_DIR, "targeted_multihabitat_scenario_summary.csv")
overall_path <- file.path(OUT_DIR, "targeted_multihabitat_overall_summary.csv")
beta_path <- file.path(OUT_DIR, "targeted_multihabitat_by_beta_summary.csv")
design_path <- file.path(OUT_DIR, "targeted_multihabitat_by_design_summary.csv")
diag_path <- file.path(OUT_DIR, "targeted_multihabitat_stan_diagnostics.csv")

readr::write_csv(results_df, results_path)
readr::write_csv(scenario_summary, scenario_path)
readr::write_csv(overall_summary, overall_path)
readr::write_csv(beta_summary, beta_path)
readr::write_csv(design_summary, design_path)
readr::write_csv(diag_df, diag_path)


# ============================================================
# 14. PRINT OUTPUT
# ============================================================

cat("\n==============================\n")
cat("OVERALL METHOD PERFORMANCE\n")
cat("==============================\n")
print(overall_summary, n = Inf)

cat("\n==============================\n")
cat("PERFORMANCE BY TRUE BETA\n")
cat("==============================\n")
print(beta_summary, n = Inf)

cat("\n==============================\n")
cat("PERFORMANCE BY H / SAMPLE SIZE DESIGN\n")
cat("==============================\n")
print(design_summary, n = Inf)

cat("\n==============================\n")
cat("SCENARIO-SPECIFIC PERFORMANCE\n")
cat("==============================\n")
print(scenario_summary, n = Inf)

cat("\n==============================\n")
cat("STAN DIAGNOSTICS\n")
cat("==============================\n")
print(as.data.frame(diag_df))

cat("\nSaved files:\n")
cat("  ", results_path, "\n")
cat("  ", scenario_path, "\n")
cat("  ", overall_path, "\n")
cat("  ", beta_path, "\n")
cat("  ", design_path, "\n")
cat("  ", diag_path, "\n")

cat("\nInterpretation notes:\n")
cat("  This is a targeted extension, not a replacement for the main 20-method benchmark.\n")
cat("  The focus is the multi-habitat regime where EIV and habitat clustering coexist.\n")
cat("  The Bayesian model explicitly estimates latent x_star and habitat random effects.\n")
cat("  Uniform priors are implemented as broad bounded parameters in Stan.\n")
cat("  Check Rhat near 1.00, adequate ESS, and low divergent transitions.\n")
cat("  After this pilot works, increase N_REPLICATES to 20, 50, or 100.\n")
