# ============================================================================
# POWER-LAW SCALING UNDER MEASUREMENT ERROR - BATCH PROCESSING FOR ALL SCENARIOS
# WITH TENSORFLOW/KERAS3 SUPPORT
# ============================================================================
# This version runs ALL scenarios automatically:
# - Different measurement error combinations: (20,40), (30,30), (40,20) as percentages
# - Different sample sizes: n=36 and n=72 per habitat
# - Different number of habitats: 6 and 12
# ============================================================================

rm(list = ls())

# ============================================================================
# LOAD REQUIRED PACKAGES
# ============================================================================

# List of packages needed
packages_needed <- c(
  "tidyverse", "parallel", "doParallel", "lme4", "smatr", "nlme",
  "lmodel2", "deming", "MethComp", "mblm", "odr", "simex", "brms", 
  "rstan", "posterior", "bayesplot", "tidybayes", "caret", "xgboost", 
  "randomForest", "glmnet", "MASS", "future", "furrr", "tictoc", 
  "progressr", "gridExtra", "cowplot", "mcr", "nnet", "reticulate", 
  "tensorflow", "keras3"
)

# ONLY LOAD packages that are already installed
for (pkg in packages_needed) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓ Loaded:", pkg, "\n")
  } else {
    cat("⚠ Not installed:", pkg, "- some features may be disabled\n")
  }
}

# ============================================================================
# CHECK KERAS/TENSORFLOW AVAILABILITY - FIXED VERSION
# ============================================================================

has_keras <- FALSE
has_tensorflow <- FALSE

# Check if tensorflow is working
if (require(tensorflow, quietly = TRUE)) {
  tryCatch({
    # Test if TensorFlow works
    tf <- tensorflow::tf
    test_result <- tf$constant("Testing TensorFlow")
    has_tensorflow <- TRUE
    
    if (require(keras3, quietly = TRUE)) {
      has_keras <- TRUE
      cat("\n✅ TensorFlow/Keras3 is WORKING!\n")
      cat("   TensorFlow version:", tf$`__version__`, "\n")
      
      # Try different ways to get keras version
      keras_ver <- tryCatch({
        keras3::keras_version()
      }, error = function(e) {
        tryCatch({
          keras3::version()
        }, error = function(e) {
          "unknown (version function not found)"
        })
      })
      cat("   Keras version:", keras_ver, "\n")
      
      # Check for GPU
      gpus <- tryCatch({
        tf$config$list_physical_devices("GPU")
      }, error = function(e) {
        list()
      })
      if (length(gpus) > 0) {
        cat("   GPU available\n")
      } else {
        cat("   Using CPU\n")
      }
    }
  }, error = function(e) {
    cat("\n❌ TensorFlow test failed:", e$message, "\n")
    cat("   Keras methods will be disabled\n")
    has_keras <- FALSE
    has_tensorflow <- FALSE
  })
} else {
  cat("\n⚠ TensorFlow not loaded properly\n")
}

# ============================================================================
# SET UP PARALLEL COMPUTING
# ============================================================================

n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = min(n_cores, 4))
options(mc.cores = min(n_cores, 4))
options(brms.backend = "cmdstanr")
Sys.setenv(OMP_NUM_THREADS = 1)

# ============================================================================
# SIMULATION PARAMETERS (COMMON ACROSS ALL SCENARIOS)
# ============================================================================

# Simulation parameters
n_sim <- 100  # Number of simulations per scenario (reduced for testing, use 300 for final)

# True parameter values
beta_values <- c(0.6, 1.0, 1.4)  # True power-law exponents to test
alpha_true <- 1.5  # True intercept (scaling factor)

# Random effects parameters
sd_beta <- 0.40  # SD for random slopes
sd_alpha <- 0.40  # SD for random intercepts
cor_alpha_beta <- -0.3  # Correlation between intercepts and slopes

# X variable parameters
x_min <- 0.5
x_max <- 4
x_means_habitat_6 <- c(1.0, 2.5, 4.0, 5.5, 7.0, 8.5)  # Habitat-specific X means for 6 habitats
x_means_habitat_12 <- c(1.0, 1.8, 2.6, 3.4, 4.2, 5.0, 5.8, 6.6, 7.4, 8.2, 9.0, 9.8)  # For 12 habitats
x_sd_habitat <- 1.2  # Spread of X within each habitat

# Process error
sigma_proc <- 0.30  # Process error
heteroscedastic <- FALSE  # Whether errors are heteroscedastic

# Bayesian parameters
brms_chains <- 2  # Reduced for testing (use 4 for final)
brms_iter <- 1000  # Reduced for testing (use 4000 for final)
brms_warmup <- 500  # Reduced for testing (use 2000 for final)
brms_adapt_delta <- 0.99  # Target acceptance rate
brms_treedepth <- 15  # Max tree depth

# Deep learning parameters
keras_epochs <- 50
keras_batch_size <- 32
keras_validation_split <- 0.2
keras_patience <- 10

# ============================================================================
# DEFINE ALL SCENARIOS TO RUN
# ============================================================================

# Base directory
base_dir <- "D:/LMEM_2026/Mevan"

# Define all scenario combinations
scenarios <- expand.grid(
  sigma_x_me = c(0.20, 0.30, 0.40),  # Measurement error in X
  sigma_y_me = c(0.20, 0.30, 0.40),  # Measurement error in Y
  n_per_habitat = c(36, 72),          # Observations per habitat
  n_habitats = c(6, 12),              # Number of habitats
  stringsAsFactors = FALSE
)

# Filter to only the specific combinations you want
# (20,40), (30,30), (40,20) for error combinations
scenarios <- scenarios[
  (scenarios$sigma_x_me == 0.20 & scenarios$sigma_y_me == 0.40) |
  (scenarios$sigma_x_me == 0.30 & scenarios$sigma_y_me == 0.30) |
  (scenarios$sigma_x_me == 0.40 & scenarios$sigma_y_me == 0.20), 
]

# Sort scenarios for logical order
scenarios <- scenarios[order(scenarios$n_habitats, 
                             scenarios$n_per_habitat, 
                             scenarios$sigma_x_me), ]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

count_divergences <- function(brmsfit_obj) {
  if (is.null(brmsfit_obj)) return(NA_real_)
  tryCatch({
    sp <- rstan::get_sampler_params(brmsfit_obj$fit, inc_warmup = FALSE)
    sum(vapply(sp, function(mat) sum(mat[, "divergent__"]), numeric(1)))
  }, error = function(e) NA_real_)
}

calc_coverage <- function(estimate, se, true_value, df = Inf) {
  if (is.na(estimate) || is.na(se)) return(NA_real_)
  ci <- if (df == Inf) {
    c(estimate - 1.96 * se, estimate + 1.96 * se)
  } else {
    estimate + c(-1, 1) * qt(0.975, df) * se
  }
  as.numeric(true_value >= ci[1] & true_value <= ci[2])
}

create_output_dir <- function(sigma_x_me, sigma_y_me, n_per_habitat, n_habitats) {
  # Convert to percentages for folder name
  x_percent <- sigma_x_me * 100
  y_percent <- sigma_y_me * 100
  
  dir_name <- sprintf("habitat_visualization_output_R_%.0f_%.0f_n%d_%d",
                      x_percent, y_percent, n_per_habitat, n_habitats)
  
  full_path <- file.path(base_dir, dir_name)
  
  # Create directory if it doesn't exist
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    cat("  Created directory:", full_path, "\n")
  }
  
  return(full_path)
}

# ============================================================================
# MAIN SIMULATION FUNCTION
# ============================================================================

run_sim <- function(beta_true, sim_id, sigma_x_me, sigma_y_me, 
                    n_per_habitat, n_habitats, multi_habitat = TRUE) {
  
  set.seed(20000 + sim_id * 1000 + beta_true * 100 + 
           sigma_x_me * 1000 + sigma_y_me * 100 + 
           n_per_habitat * 10 + n_habitats)
  
  # ==========================================================================
  # DATA GENERATION WITH RANDOM SLOPES AND HABITAT-SPECIFIC X MEANS
  # ==========================================================================
  
  # Select appropriate X means based on number of habitats
  if (n_habitats == 6) {
    x_means_habitat <- c(1.0, 2.5, 4.0, 5.5, 7.0, 8.5)
  } else {
    x_means_habitat <- c(1.0, 1.8, 2.6, 3.4, 4.2, 5.0, 5.8, 6.6, 7.4, 8.2, 9.0, 9.8)
  }
  
  if (multi_habitat) {
    habitat <- rep(1:n_habitats, each = n_per_habitat)
    N <- length(habitat)
    
    # Generate correlated random effects (intercepts and slopes)
    Sigma <- matrix(c(sd_alpha^2, cor_alpha_beta * sd_alpha * sd_beta,
                      cor_alpha_beta * sd_alpha * sd_beta, sd_beta^2), 2, 2)
    
    # Generate random effects
    if (requireNamespace("MASS", quietly = TRUE)) {
      random_effects <- MASS::mvrnorm(n_habitats, mu = c(0, 0), Sigma = Sigma)
    } else {
      # Fallback if MASS not available
      random_effects <- cbind(
        rnorm(n_habitats, 0, sd_alpha),
        rnorm(n_habitats, 0, sd_beta)
      )
    }
    
    # Habitat-specific intercepts and slopes
    beta_h <- beta_true + random_effects[, 2]
    alpha_h <- pmax(alpha_true + random_effects[, 1], 0.1)  # Ensure positive
    
    # Generate X values with habitat-specific means
    log_x_true <- vector(length = N)
    for (h in 1:n_habitats) {
      idx <- which(habitat == h)
      # Generate log(X) with habitat-specific mean
      log_mean <- log(x_means_habitat[h])
      log_sd <- log(x_sd_habitat + 1) / 2  # Adjust for log scale
      log_x_true[idx] <- rnorm(n_per_habitat, mean = log_mean, sd = log_sd)
    }
    
    # Generate Y values with habitat-specific parameters
    log_y_true <- vector(length = N)
    for (i in 1:N) {
      h <- habitat[i]
      log_y_true[i] <- log(alpha_h[h]) + beta_h[h] * log_x_true[i] + rnorm(1, 0, sigma_proc)
    }
    
  } else {
    # Single habitat case (fallback)
    habitat <- rep(1, n_per_habitat)
    N <- length(habitat)
    beta_h <- beta_true
    alpha_h <- alpha_true
    log_x_true <- runif(N, log(x_min), log(x_max))
    log_y_true <- log(alpha_h) + beta_h * log_x_true + rnorm(N, 0, sigma_proc)
  }
  
  # Add measurement error
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
    habitat = factor(habitat), 
    log_x = log_x_obs, 
    log_y = log_y_obs,
    se_x = se_x, 
    se_y = se_y, 
    x_true = exp(log_x_true), 
    y_true = exp(log_y_true),
    x_obs = exp(log_x_obs), 
    y_obs = exp(log_y_obs),
    alpha_true_habitat = alpha_h[habitat],
    beta_true_habitat = beta_h[habitat]
  )
  
  # Initialize results
  results <- list(
    sim_id = sim_id, 
    beta_true = beta_true, 
    n = N,
    sigma_x_me = sigma_x_me,
    sigma_y_me = sigma_y_me,
    n_per_habitat = n_per_habitat,
    n_habitats = n_habitats,
    
    # CLASSICAL METHODS
    beta_ols = NA_real_, se_ols = NA_real_, cover_ols = NA_real_,
    beta_sma = NA_real_,
    beta_majoraxis = NA_real_,
    beta_rma = NA_real_,
    beta_odr = NA_real_,
    
    # DEMING FAMILY
    beta_deming_std = NA_real_, se_deming_std = NA_real_, cover_deming_std = NA_real_,
    beta_deming_wtd = NA_real_,
    beta_deming_mcr = NA_real_, 
    lambda_estimated = NA_real_,
    beta_pbablok = NA_real_,
    
    # NON-PARAMETRIC
    beta_theilsen = NA_real_,
    beta_siegel = NA_real_,
    
    # SIMEX
    beta_simex = NA_real_, se_simex = NA_real_,
    
    # MIXED MODELS
    beta_lmer = NA_real_, se_lmer = NA_real_, cover_lmer = NA_real_,
    beta_lmer_slopes = NA_real_, se_lmer_slopes = NA_real_,
    beta_nlme = NA_real_,
    
    # BAYESIAN
    beta_brms_std = NA_real_, se_brms_std = NA_real_, cover_brms_std = NA_real_,
    beta_brms_me = NA_real_, se_brms_me = NA_real_, cover_brms_me = NA_real_,
    beta_brms_robust = NA_real_, se_brms_robust = NA_real_,
    beta_brms_horseshoe = NA_real_,
    
    # ML METHODS
    beta_xgboost = NA_real_,
    beta_xgboost_corrected = NA_real_,
    beta_rf = NA_real_,
    beta_glmnet = NA_real_,
    beta_nnet = NA_real_,
    beta_keras = NA_real_,
    beta_keras_denoised = NA_real_,
    
    # ENSEMBLE
    beta_ensemble = NA_real_,
    
    # Diagnostics
    rhat = NA_real_, ess_bulk = NA_real_, divergent = NA_real_,
    time_total = NA_real_,
    
    sd_beta_actual = sd(beta_h),
    sd_alpha_actual = sd(alpha_h),
    cor_alpha_beta_actual = if(length(alpha_h) > 1) cor(alpha_h, beta_h) else NA_real_
  )
  
  sim_start <- Sys.time()
  
  # ==========================================================================
  # SECTION 1: CLASSICAL METHODS
  # ==========================================================================
  
  # 1.1 OLS
  tryCatch({
    m_ols <- lm(log_y ~ log_x, data = df)
    results$beta_ols <- unname(coef(m_ols)[2])
    results$se_ols <- summary(m_ols)$coefficients[2, 2]
    ci <- confint(m_ols)[2,]
    results$cover_ols <- as.numeric(beta_true >= ci[1] & beta_true <= ci[2])
  }, error = function(e) {})
  
  # 1.2 SMA
  tryCatch({
    m_sma <- sma(log_y ~ log_x, data = df)
    if (!is.null(m_sma$coef)) results$beta_sma <- as.numeric(m_sma$coef[[1]]["slope"])
  }, error = function(e) {})
  
  # 1.3 Major Axis and RMA
  tryCatch({
    ma_reg <- suppressMessages(lmodel2(log_y ~ log_x, data = df, 
                                       range.y = "interval", range.x = "interval", nperm = 0))
    if (!is.null(ma_reg$regression.results)) {
      results$beta_majoraxis <- ma_reg$regression.results$Slope[2]
      results$beta_rma <- ma_reg$regression.results$Slope[3]
    }
  }, error = function(e) {})
  
  # 1.4 Orthogonal Distance Regression
  tryCatch({
    odr_fit <- odr(log_y ~ log_x, data = df, delta = se_x, eps = se_y)
    results$beta_odr <- coef(odr_fit)[2]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 2: DEMING FAMILY
  # ==========================================================================
  
  # 2.1 Standard Deming
  tryCatch({
    dem_fit <- deming(log_y ~ log_x, data = df, xstd = se_x, ystd = se_y)
    results$beta_deming_std <- coef(dem_fit)[2]
    results$se_deming_std <- sqrt(diag(vcov(dem_fit)))[2]
    ci_lower <- results$beta_deming_std - 1.96 * results$se_deming_std
    ci_upper <- results$beta_deming_std + 1.96 * results$se_deming_std
    results$cover_deming_std <- as.numeric(beta_true >= ci_lower & beta_true <= ci_upper)
  }, error = function(e) {})
  
  # 2.2 Weighted Deming
  tryCatch({
    wts <- 1 / (se_x^2 + se_y^2)
    dem_wtd <- deming(log_y ~ log_x, data = df, xstd = se_x, ystd = se_y, weights = wts)
    results$beta_deming_wtd <- coef(dem_wtd)[2]
  }, error = function(e) {})
  
  # 2.3 Deming with mcr package
  tryCatch({
    dem_mcr <- mcreg(df$log_x, df$log_y, method.reg = "Deming", 
                     error.ratio = (sigma_y_me^2)/(sigma_x_me^2))
    results$beta_deming_mcr <- dem_mcr@para[2, 1]
    results$lambda_estimated <- (sigma_y_me^2)/(sigma_x_me^2)
  }, error = function(e) {})
  
  # 2.4 Passing-Bablok
  tryCatch({
    pb_reg <- mcreg(df$log_x, df$log_y, method.reg = "PaBa")
    results$beta_pbablok <- pb_reg@para[2, 1]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 3: NON-PARAMETRIC
  # ==========================================================================
  
  # 3.1 Theil-Sen
  tryCatch({
    ts_reg <- mblm(log_y ~ log_x, df, repeated = TRUE)
    results$beta_theilsen <- coef(ts_reg)[2]
  }, error = function(e) {})
  
  # 3.2 Siegel
  tryCatch({
    siegel_reg <- mblm(log_y ~ log_x, df, repeated = FALSE)
    results$beta_siegel <- coef(siegel_reg)[2]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 4: SIMEX
  # ==========================================================================
  
  tryCatch({
    naive_fit <- lm(log_y ~ log_x, data = df)
    simex_fit <- simex(naive_fit, 
                       SIMEXvariable = "log_x",
                       measurement.error = sigma_x_me,
                       lambda = seq(0, 2, 0.2),
                       extrapolation = "quadratic",
                       B = 50)
    results$beta_simex <- coef(simex_fit)[2]
    results$se_simex <- sqrt(diag(vcov(simex_fit)))[2]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 5: MIXED MODELS
  # ==========================================================================
  
  if (multi_habitat) {
    # 5.1 Linear Mixed Model (random intercepts only)
    tryCatch({
      m_lmer <- lmer(log_y ~ log_x + (1 | habitat), data = df, 
                     REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      results$beta_lmer <- unname(fixef(m_lmer)[2])
      results$se_lmer <- sqrt(vcov(m_lmer)[2, 2])
      ci_lower <- results$beta_lmer - 1.96 * results$se_lmer
      ci_upper <- results$beta_lmer + 1.96 * results$se_lmer
      results$cover_lmer <- as.numeric(beta_true >= ci_lower & beta_true <= ci_upper)
    }, error = function(e) {})
    
    # 5.2 Random slopes model
    tryCatch({
      m_lmer_slopes <- lmer(log_y ~ log_x + (log_x | habitat), data = df, 
                             REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      results$beta_lmer_slopes <- unname(fixef(m_lmer_slopes)[2])
      results$se_lmer_slopes <- sqrt(vcov(m_lmer_slopes)[2, 2])
    }, error = function(e) {})
  }
  
  # 5.3 Nonlinear Model
  tryCatch({
    m_nls <- nls(y_obs ~ a * x_obs^b, data = df,
                 start = c(a = alpha_true, b = beta_true),
                 control = nls.control(maxiter = 100))
    results$beta_nlme <- coef(m_nls)["b"]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 6: BAYESIAN METHODS
  # ==========================================================================
  
  # 6.1 BRMS Standard
  tryCatch({
    m_brms_std <- brm(
      bf(log_y ~ log_x + (1 | habitat)),
      data = df,
      prior = c(
        prior(normal(0, 1), class = "b"),
        prior(normal(1, 0.5), class = "Intercept"),
        prior(exponential(5), class = "sd")
      ),
      chains = brms_chains,
      iter = brms_iter,
      warmup = brms_warmup,
      control = list(adapt_delta = brms_adapt_delta, max_treedepth = brms_treedepth),
      refresh = 0, 
      silent = 2
    )
    results$beta_brms_std <- fixef(m_brms_std)["log_x", "Estimate"]
    results$se_brms_std <- fixef(m_brms_std)["log_x", "Est.Error"]
    post_int <- posterior_interval(m_brms_std, prob = 0.95)
    results$cover_brms_std <- as.numeric(beta_true >= post_int["log_x", 1] & 
                                         beta_true <= post_int["log_x", 2])
    results$rhat <- max(rhat(m_brms_std))
    results$ess_bulk <- min(ess_bulk(m_brms_std))
    results$divergent <- count_divergences(m_brms_std)
  }, error = function(e) {})
  
  # 6.2 BRMS with Measurement Error
  tryCatch({
    m_brms_me <- brm(
      bf(log_y | se(se_y) ~ 1 + me(log_x, se_x) + (1 | habitat)),
      data = df,
      prior = c(
        prior(normal(0, 2), class = "b"),
        prior(normal(1, 1), class = "Intercept"),
        prior(exponential(2), class = "sd"),
        prior(gamma(4, 2), class = "sderr")
      ),
      chains = brms_chains,
      iter = brms_iter,
      warmup = brms_warmup,
      control = list(adapt_delta = 0.999, max_treedepth = 18),
      refresh = 0, 
      silent = 2
    )
    me_name <- grep("^me", rownames(fixef(m_brms_me)), value = TRUE)[1]
    if (!is.na(me_name)) {
      results$beta_brms_me <- fixef(m_brms_me)[me_name, "Estimate"]
      results$se_brms_me <- fixef(m_brms_me)[me_name, "Est.Error"]
      post_int <- posterior_interval(m_brms_me, prob = 0.95)
      if (me_name %in% rownames(post_int)) {
        results$cover_brms_me <- as.numeric(beta_true >= post_int[me_name, 1] & 
                                            beta_true <= post_int[me_name, 2])
      }
    }
  }, error = function(e) {})
  
  # 6.3 BRMS Robust
  tryCatch({
    m_brms_robust <- brm(
      bf(log_y ~ log_x + (1 | habitat)),
      data = df, 
      family = student(),
      prior = c(
        prior(normal(0, 1), class = "b"),
        prior(normal(1, 0.5), class = "Intercept"),
        prior(exponential(5), class = "sd"),
        prior(gamma(2, 0.1), class = "nu")
      ),
      chains = brms_chains,
      iter = brms_iter,
      warmup = brms_warmup,
      control = list(adapt_delta = brms_adapt_delta, max_treedepth = brms_treedepth),
      refresh = 0, 
      silent = 2
    )
    results$beta_brms_robust <- fixef(m_brms_robust)["log_x", "Estimate"]
    results$se_brms_robust <- fixef(m_brms_robust)["log_x", "Est.Error"]
  }, error = function(e) {})
  
  # 6.4 BRMS Horseshoe
  tryCatch({
    m_brms_hs <- brm(
      bf(log_y ~ log_x + (1 | habitat)),
      data = df,
      prior = c(
        prior(horseshoe(3), class = "b"),
        prior(normal(1, 0.5), class = "Intercept"),
        prior(exponential(5), class = "sd")
      ),
      chains = brms_chains,
      iter = brms_iter,
      warmup = brms_warmup,
      control = list(adapt_delta = brms_adapt_delta, max_treedepth = brms_treedepth),
      refresh = 0, 
      silent = 2
    )
    results$beta_brms_horseshoe <- fixef(m_brms_hs)["log_x", "Estimate"]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 7: MACHINE LEARNING METHODS
  # ==========================================================================
  
  # 7.1 XGBoost
  tryCatch({
    X <- as.matrix(df[, c("log_x", "se_x", "se_y")])
    y <- df$log_y
    dtrain <- xgb.DMatrix(X, label = y)
    params <- list(objective = "reg:squarederror", eta = 0.1, max_depth = 3)
    xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 50, verbose = 0)
    importance <- xgb.importance(model = xgb_model)
    log_x_imp <- importance$Gain[importance$Feature == "log_x"]
    if (length(log_x_imp) > 0) results$beta_xgboost <- log_x_imp
  }, error = function(e) {})
  
  # 7.2 XGBoost Corrected
  tryCatch({
    X_corr <- as.matrix(df[, c("log_x", "se_x", "se_y", "I(log_x^2)")])
    dtrain_corr <- xgb.DMatrix(X_corr, label = y)
    xgb_corr <- xgb.train(params = params, data = dtrain_corr, nrounds = 50, verbose = 0)
    imp_corr <- xgb.importance(model = xgb_corr)
    log_x_corr <- imp_corr$Gain[imp_corr$Feature == "log_x"]
    if (length(log_x_corr) > 0) results$beta_xgboost_corrected <- log_x_corr
  }, error = function(e) {})
  
  # 7.3 Random Forest
  tryCatch({
    rf_data <- df[, c("log_y", "log_x", "se_x", "se_y")]
    rf_model <- randomForest(log_y ~ ., data = rf_data, ntree = 100)
    pd <- partialPlot(rf_model, rf_data, "log_x", plot = FALSE)
    if (!is.null(pd) && length(pd$x) > 1) {
      results$beta_rf <- coef(lm(pd$y ~ pd$x))[2]
    }
  }, error = function(e) {})
  
  # 7.4 glmnet
  tryCatch({
    x_mat <- model.matrix(~ log_x + I(log_x^2) + habitat, data = df)[, -1]
    cv_fit <- cv.glmnet(x_mat, df$log_y, alpha = 0)
    results$beta_glmnet <- coef(cv_fit, s = "lambda.min")["log_x", 1]
  }, error = function(e) {})
  
  # 7.5 Neural Network (nnet)
  tryCatch({
    nn_fit <- nnet(log_y ~ log_x + se_x + se_y, data = df, size = 3, 
                   linout = TRUE, trace = FALSE, maxit = 100)
    results$beta_nnet <- coef(nn_fit)[2]
  }, error = function(e) {})
  
  # ==========================================================================
  # SECTION 8: KERAS3 NEURAL NETWORK METHODS
  # ==========================================================================
  
  if (has_keras) {
    
    # 8.1 Simple Neural Network with Keras3
    tryCatch({
      # Prepare data
      x_train <- as.matrix(df[, c("log_x", "se_x", "se_y")])
      y_train <- as.matrix(df$log_y)
      
      # Build a simple neural network
      model_keras <- keras_model_sequential(input_shape = ncol(x_train)) %>%
        layer_dense(units = 16, activation = "relu") %>%
        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 8, activation = "relu") %>%
        layer_dense(units = 1)
      
      # Compile model
      model_keras %>% compile(
        optimizer = optimizer_adam(learning_rate = 0.001),
        loss = "mse",
        metrics = c("mae")
      )
      
      # Early stopping callback
      early_stop <- callback_early_stopping(
        monitor = "val_loss",
        patience = keras_patience,
        restore_best_weights = TRUE
      )
      
      # Train model
      history <- model_keras %>% fit(
        x_train, y_train,
        epochs = keras_epochs,
        batch_size = keras_batch_size,
        validation_split = keras_validation_split,
        verbose = 0,
        callbacks = list(early_stop)
      )
      
      # Extract feature importance
      weights <- get_weights(model_keras)
      if (length(weights) >= 2) {
        first_layer_weights <- weights[[1]]
        results$beta_keras <- mean(abs(first_layer_weights[1, ]))
      }
      
    }, error = function(e) {})
    
    # 8.2 Denoising Autoencoder with Keras3
    tryCatch({
      # Prepare data
      X_ae <- as.matrix(cbind(df$log_x, df$log_y, df$se_x, df$se_y))
      
      # Build autoencoder
      input_dim <- ncol(X_ae)
      encoding_dim <- 2
      
      # Encoder
      input_layer <- layer_input(shape = c(input_dim))
      encoded <- input_layer %>%
        layer_dense(units = 8, activation = "relu") %>%
        layer_dense(units = encoding_dim, activation = "linear")
      
      # Decoder
      decoded <- encoded %>%
        layer_dense(units = 8, activation = "relu") %>%
        layer_dense(units = input_dim, activation = "linear")
      
      # Autoencoder model
      autoencoder <- keras_model(inputs = input_layer, outputs = decoded)
      
      # Compile
      autoencoder %>% compile(
        optimizer = optimizer_adam(learning_rate = 0.001),
        loss = "mse"
      )
      
      # Early stopping
      early_stop_ae <- callback_early_stopping(
        monitor = "val_loss",
        patience = keras_patience,
        restore_best_weights = TRUE
      )
      
      # Train autoencoder
      autoencoder %>% fit(
        X_ae, X_ae,
        epochs = keras_epochs,
        batch_size = keras_batch_size,
        validation_split = keras_validation_split,
        verbose = 0,
        callbacks = list(early_stop_ae)
      )
      
      # Get denoised data
      X_denoised <- predict(autoencoder, X_ae)
      
      # Create denoised dataframe
      df_denoised <- df
      df_denoised$log_x_denoised <- X_denoised[, 1]
      df_denoised$log_y_denoised <- X_denoised[, 2]
      
      # Fit linear model on denoised data
      m_denoised <- lm(log_y_denoised ~ log_x_denoised, data = df_denoised)
      results$beta_keras_denoised <- coef(m_denoised)[2]
      
    }, error = function(e) {})
  }
  
  # ==========================================================================
  # SECTION 9: ENSEMBLE
  # ==========================================================================
  
  tryCatch({
    available <- c(
      results$beta_ols, 
      results$beta_deming_std, 
      results$beta_theilsen, 
      results$beta_majoraxis,
      results$beta_brms_std,
      results$beta_rf,
      results$beta_keras,
      results$beta_keras_denoised
    )
    available <- available[!is.na(available) & !is.infinite(available) & abs(available) < 10]
    if (length(available) >= 3) {
      results$beta_ensemble <- mean(available, trim = 0.2)
    }
  }, error = function(e) {})
  
  # Record time
  results$time_total <- as.numeric(difftime(Sys.time(), sim_start, units = "secs"))
  
  return(as_tibble(results))
}

# ============================================================================
# FUNCTION TO RUN ALL SIMULATIONS FOR A SINGLE SCENARIO
# ============================================================================

run_scenario <- function(sigma_x_me, sigma_y_me, n_per_habitat, n_habitats) {
  
  # Create output directory
  out_dir <- create_output_dir(sigma_x_me, sigma_y_me, n_per_habitat, n_habitats)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("RUNNING SCENARIO:\n")
  cat(sprintf("  σx = %.2f, σy = %.2f\n", sigma_x_me, sigma_y_me))
  cat(sprintf("  n per habitat = %d, n habitats = %d\n", n_per_habitat, n_habitats))
  cat(sprintf("  Output directory: %s\n", out_dir))
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # Create parameter grid for this scenario
  sim_params <- expand.grid(
    beta_true = beta_values, 
    sim_id = 1:n_sim,
    sigma_x_me = sigma_x_me,
    sigma_y_me = sigma_y_me,
    n_per_habitat = n_per_habitat,
    n_habitats = n_habitats
  )
  
  cat(sprintf("Starting %d simulations (%d β values × %d reps)...\n", 
              nrow(sim_params), length(beta_values), n_sim))
  
  # Run simulations with progress bar
  pb <- txtProgressBar(min = 0, max = nrow(sim_params), style = 3)
  results_list <- list()
  
  for (i in 1:nrow(sim_params)) {
    results_list[[i]] <- run_sim(
      beta_true = sim_params$beta_true[i],
      sim_id = sim_params$sim_id[i],
      sigma_x_me = sigma_x_me,
      sigma_y_me = sigma_y_me,
      n_per_habitat = n_per_habitat,
      n_habitats = n_habitats
    )
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combine results
  all_results <- bind_rows(results_list)
  
  cat(sprintf("\n✅ Simulations complete for this scenario! Total runs: %d\n", nrow(all_results)))
  
  # ==========================================================================
  # POST-PROCESSING
  # ==========================================================================
  
  cat("\n📊 Post-processing results...\n")
  
  # Filter out extreme values
  results_clean <- all_results %>%
    mutate(across(starts_with("beta_"), ~ifelse(abs(.x) > 10, NA, .x)))
  
  # Find working methods
  method_success <- results_clean %>%
    summarise(across(starts_with("beta_"), ~sum(!is.na(.x)), .names = "{.col}")) %>%
    pivot_longer(everything(), names_to = "method", values_to = "n_success") %>%
    mutate(method = str_remove(method, "beta_")) %>%
    filter(n_success > 0, method != "true") %>%
    arrange(desc(n_success))
  
  cat("\n📊 Working methods:\n")
  print(method_success, n = Inf)
  
  # ==========================================================================
  # CREATE FIGURES
  # ==========================================================================
  
  cat("\n📈 Creating figures...\n")
  
  if (nrow(method_success) > 0) {
    
    # Remove nnet from plotting methods
    plot_methods <- setdiff(method_success$method, "nnet")
    plot_methods <- plot_methods[method_success$n_success[method_success$method %in% plot_methods] > 10]
    
    cat("Plotting methods:", paste(plot_methods, collapse = ", "), "\n")
    
    # FIGURE 1: Boxplots
    plot_data <- results_clean %>%
      dplyr::select(beta_true, all_of(paste0("beta_", plot_methods))) %>%
      tidyr::pivot_longer(cols = -beta_true, names_to = "method", values_to = "estimate") %>%
      mutate(
        method = str_remove(method, "beta_"),
        beta_true = factor(beta_true)
      ) %>%
      filter(!is.na(estimate))
    
    p1 <- ggplot(plot_data, aes(x = reorder(method, estimate, median), y = estimate, fill = method)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_hline(data = distinct(plot_data, beta_true),
                 aes(yintercept = as.numeric(as.character(beta_true))), 
                 linetype = "dashed", color = "red", linewidth = 0.8) +
      facet_wrap(~beta_true, ncol = 3, 
                 labeller = labeller(beta_true = function(x) paste("True β =", x))) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "lightgray")
      ) +
      labs(
        title = sprintf("Method Comparison: σx=%.2f, σy=%.2f, n=%d, habitats=%d",
                        sigma_x_me, sigma_y_me, n_per_habitat, n_habitats),
        x = "Method", y = "Estimated Exponent (β̂)"
      )
    
    ggsave(file.path(out_dir, "Figure1_boxplots.png"), p1, width = 14, height = 8, dpi = 300)
    
    # FIGURE 2: Bias Bar Plot
    bias_data <- results_clean %>%
      dplyr::select(beta_true, all_of(paste0("beta_", plot_methods))) %>%
      tidyr::pivot_longer(cols = -beta_true, names_to = "method", values_to = "estimate") %>%
      mutate(
        method = str_remove(method, "beta_"),
        bias = estimate - beta_true
      ) %>%
      group_by(beta_true, method) %>%
      summarise(
        mean_bias = mean(bias, na.rm = TRUE),
        se_bias = sd(bias, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    p2 <- ggplot(bias_data, aes(x = method, y = mean_bias, fill = factor(beta_true))) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = mean_bias - se_bias, ymax = mean_bias + se_bias),
                    position = position_dodge(width = 0.8), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      ) +
      labs(
        title = sprintf("Method Bias: σx=%.2f, σy=%.2f, n=%d, habitats=%d",
                        sigma_x_me, sigma_y_me, n_per_habitat, n_habitats),
        x = "Method", y = "Mean Bias (β̂ - β)",
        fill = "True β"
      )
    
    ggsave(file.path(out_dir, "Figure2_bias.png"), p2, width = 14, height = 8, dpi = 300)
    
    # FIGURE 3: RMSE Heatmap
    rmse_data <- results_clean %>%
      dplyr::select(beta_true, all_of(paste0("beta_", plot_methods))) %>%
      tidyr::pivot_longer(cols = -beta_true, names_to = "method", values_to = "estimate") %>%
      mutate(
        method = str_remove(method, "beta_"),
        se = (estimate - beta_true)^2
      ) %>%
      group_by(beta_true, method) %>%
      summarise(
        rmse = sqrt(mean(se, na.rm = TRUE)),
        .groups = "drop"
      )
    
    p3 <- ggplot(rmse_data, aes(x = method, y = factor(beta_true), fill = rmse)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(rmse, 3)), size = 3.5, color = "black") +
      scale_fill_gradient(low = "white", high = "red", name = "RMSE") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      labs(
        title = sprintf("RMSE Heatmap: σx=%.2f, σy=%.2f, n=%d, habitats=%d",
                        sigma_x_me, sigma_y_me, n_per_habitat, n_habitats),
        x = "Method", y = "True β"
      )
    
    ggsave(file.path(out_dir, "Figure3_rmse_heatmap.png"), p3, width = 14, height = 6, dpi = 300)
    
    # FIGURE 4: Method Ranking
    ranking_data <- rmse_data %>%
      group_by(method) %>%
      summarise(avg_rmse = mean(rmse)) %>%
      arrange(avg_rmse) %>%
      mutate(rank = row_number())
    
    p4 <- ggplot(ranking_data, aes(x = reorder(method, -avg_rmse), y = avg_rmse, fill = avg_rmse)) +
      geom_col() +
      geom_text(aes(label = paste0("Rank ", rank, "\n", round(avg_rmse, 3))), 
                vjust = -0.5, size = 3) +
      scale_fill_gradient(low = "green", high = "red") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      labs(
        title = sprintf("Method Ranking by RMSE: σx=%.2f, σy=%.2f, n=%d, habitats=%d",
                        sigma_x_me, sigma_y_me, n_per_habitat, n_habitats),
        x = "Method", y = "Average RMSE"
      )
    
    ggsave(file.path(out_dir, "Figure4_ranking.png"), p4, width = 12, height = 6, dpi = 300)
    
    cat("✓ All figures saved\n")
  }
  
  # ==========================================================================
  # CREATE SUMMARY STATISTICS
  # ==========================================================================
  
  cat("\n📊 Creating summary statistics...\n")
  
  # Summary table by true beta
  summary_stats <- results_clean %>%
    group_by(beta_true) %>%
    summarise(
      n = n(),
      across(all_of(paste0("beta_", plot_methods)), 
             list(
               mean = ~mean(.x, na.rm = TRUE),
               sd = ~sd(.x, na.rm = TRUE),
               bias = ~mean(.x - beta_true, na.rm = TRUE),
               rmse = ~sqrt(mean((.x - beta_true)^2, na.rm = TRUE))
             ), .names = "{.col}_{.fn}")
    )
  
  # ==========================================================================
  # SAVE RESULTS
  # ==========================================================================
  
  cat("\n💾 Saving results...\n")
  
  write_csv(all_results, file.path(out_dir, "results_raw.csv"))
  write_csv(results_clean, file.path(out_dir, "results_clean.csv"))
  write_csv(method_success, file.path(out_dir, "working_methods.csv"))
  write_csv(summary_stats, file.path(out_dir, "summary_statistics.csv"))
  write_csv(ranking_data, file.path(out_dir, "method_ranking.csv"))
  
  saveRDS(list(
    raw = all_results, 
    clean = results_clean, 
    methods = method_success,
    summary = summary_stats,
    ranking = ranking_data,
    params = list(
      n_sim = n_sim, 
      sigma_x_me = sigma_x_me,
      sigma_y_me = sigma_y_me, 
      beta_values = beta_values,
      n_per_habitat = n_per_habitat,
      n_habitats = n_habitats,
      sd_beta = sd_beta,
      sd_alpha = sd_alpha,
      cor_alpha_beta = cor_alpha_beta,
      brms_chains = brms_chains,
      brms_iter = brms_iter,
      keras_epochs = keras_epochs,
      has_keras = has_keras
    )
  ), file.path(out_dir, "complete_results.rds"))
  
  # Save session info
  sink(file.path(out_dir, "session_info.txt"))
  print(sessionInfo())
  sink()
  
  cat(sprintf("\n✅ Scenario complete! Results saved to: %s\n", out_dir))
  
  # Return summary
  return(list(
    scenario = paste0("σx=", sigma_x_me, "_σy=", sigma_y_me, 
                      "_n=", n_per_habitat, "_h=", n_habitats),
    out_dir = out_dir,
    n_methods = nrow(method_success),
    top_methods = head(ranking_data$method, 5)
  ))
}

# ============================================================================
# RUN ALL SCENARIOS
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("POWER-LAW SCALING UNDER MEASUREMENT ERROR - BATCH PROCESSING\n")
cat("============================================================\n\n")

cat(sprintf("Total scenarios to run: %d\n", nrow(scenarios)))
cat("\nScenarios:\n")
for (i in 1:nrow(scenarios)) {
  cat(sprintf("  %d. σx=%.2f, σy=%.2f, n=%d, habitats=%d\n", 
              i, scenarios$sigma_x_me[i], scenarios$sigma_y_me[i], 
              scenarios$n_per_habitat[i], scenarios$n_habitats[i]))
}

cat("\nSimulation parameters:\n")
cat(sprintf("  - Simulations per scenario: %d\n", n_sim))
cat(sprintf("  - True β values: %s\n", paste(beta_values, collapse = ", ")))
cat(sprintf("  - Random slopes SD: %.2f\n", sd_beta))
cat(sprintf("  - Random intercepts SD: %.2f\n", sd_alpha))
cat(sprintf("  - Correlation (α-β): %.2f\n", cor_alpha_beta))
cat(sprintf("  - Keras3 available: %s\n", has_keras))
cat(sprintf("  - Bayesian chains: %d\n", brms_chains))
cat(sprintf("  - Bayesian iterations: %d\n", brms_iter))
cat("\n")

# Record overall start time
overall_start <- Sys.time()

# Run each scenario
scenario_summaries <- list()

for (i in 1:nrow(scenarios)) {
  cat("\n")
  cat(paste(rep("#", 70), collapse = ""), "\n")
  cat(sprintf("### SCENARIO %d OF %d\n", i, nrow(scenarios)))
  cat(paste(rep("#", 70), collapse = ""), "\n")
  
  scenario_start <- Sys.time()
  
  result <- run_scenario(
    sigma_x_me = scenarios$sigma_x_me[i],
    sigma_y_me = scenarios$sigma_y_me[i],
    n_per_habitat = scenarios$n_per_habitat[i],
    n_habitats = scenarios$n_habitats[i]
  )
  
  scenario_time <- difftime(Sys.time(), scenario_start, units = "mins")
  
  scenario_summaries[[i]] <- result
  
  cat(sprintf("\n⏱️  Scenario %d completed in %.1f minutes\n", i, scenario_time))
}

# ============================================================================
# FINAL SUMMARY REPORT
# ============================================================================

overall_time <- difftime(Sys.time(), overall_start, units = "hours")

cat("\n")
cat("============================================================\n")
cat("✅ ALL SCENARIOS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Total scenarios run: %d\n", nrow(scenarios)))
cat(sprintf("Total simulations: %d\n", nrow(scenarios) * n_sim * length(beta_values)))
cat(sprintf("Total time: %.2f hours\n", overall_time))
cat("\nScenario summaries:\n")

for (i in 1:length(scenario_summaries)) {
  cat(sprintf("\n%d. %s\n", i, scenario_summaries[[i]]$scenario))
  cat(sprintf("   Output: %s\n", scenario_summaries[[i]]$out_dir))
  cat(sprintf("   Working methods: %d\n", scenario_summaries[[i]]$n_methods))
  cat(sprintf("   Top 5 methods: %s\n", 
              paste(scenario_summaries[[i]]$top_methods, collapse = ", ")))
}

cat("\n")
cat("============================================================\n")
cat("Files saved in each scenario folder:\n")
cat("  - results_raw.csv: All simulation results\n")
cat("  - results_clean.csv: Filtered results\n")
cat("  - working_methods.csv: Methods that worked\n")
cat("  - summary_statistics.csv: Summary by β value\n")
cat("  - method_ranking.csv: Methods ranked by RMSE\n")
cat("  - complete_results.rds: Complete R object\n")
cat("  - Figure1_boxplots.png\n")
cat("  - Figure2_bias.png\n")
cat("  - Figure3_rmse_heatmap.png\n")
cat("  - Figure4_ranking.png\n")
cat("  - session_info.txt\n")
cat("============================================================\n")