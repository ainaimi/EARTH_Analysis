pacman::p_load(
  rio,
  here,
  tidyverse,
  parallel,
  future,
  future.apply,
  progressr,
  SuperLearner,
  broom,
  lmtest,
  sandwich,
  ranger,
  xgboost
)

# ==============================================================================
# SENSITIVITY ANALYSIS: Effect of Data Specifications on Linear Projection
# ==============================================================================
# This program explores how different data specifications affect the best
# linear projection test statistics for EDC effect modification.
#
# Specifications tested (16 total combinations):
# 1. Outliers: in vs out (2 options)
# 2. EDC trimming: 90th, 92.5th, 95th, 97.5th percentiles (4 options)
# 3. AFC trimming: trimmed at 30 vs not trimmed (2 options)
# ==============================================================================

set.seed(123)

# Load base data
load(here("data", "imputed_EARTH.Rdata"))
afc_base <- imputed_data %>% rename(AFCt = AFC_t)

# Load outlier IDs (flagged by 2+ methods from program 2)
outliers_id <- readRDS(here("output", "outliers_id.rds"))
outlier_rows <- outliers_id$row_id

# Define variables
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "BPA", "BP", "MP", "PP", "Hg")

the_data <- c("AFCt", "year", "month", "DOR", "eversmok",
              "age", "bmi", "races", "educ1", "smokstat", "previousIVF",
              "previousIUI", "gravid", env_vars)

impute_flags <- c("imp_sgratio_pht", "imp_smokstat", "imp_races", "imp_bmi",
                  "imp_previousIUI", "imp_age", "imp_educ1", "imp_Hg",
                  "imp_mBP", "imp_mCNP", "imp_B_PB", "imp_AFScanDate")

# ==============================================================================
# DATA PREPARATION FUNCTION
# ==============================================================================

prepare_data <- function(data,
                        remove_outliers = FALSE,
                        edc_percentile = NULL,
                        afc_trim = FALSE) {

  # Start with base data
  dat <- data[, c(the_data, impute_flags)]

  # Store row IDs before any filtering
  dat$row_id <- as.numeric(rownames(dat))

  # Remove outliers if specified
  if (remove_outliers) {
    dat <- dat %>% filter(!row_id %in% outlier_rows)
    cat("  Removed", length(outlier_rows), "outliers\n")
  }

  # Trim EDCs at specified percentile
  if (!is.null(edc_percentile)) {
    for (var in env_vars) {
      threshold <- quantile(dat[[var]], edc_percentile)
      dat[[var]][dat[[var]] > threshold] <- threshold
    }
    cat("  Trimmed EDCs at", edc_percentile * 100, "th percentile\n")
  }

  # Trim AFC at 30
  if (afc_trim) {
    dat$AFCt[dat$AFCt > 30] <- 30
    cat("  Trimmed AFC at 30\n")
  }

  # Create age binary variable
  dat$age_bin <- as.numeric(dat$age >= 35)

  # Remove row_id before returning
  dat <- dat %>% select(-row_id)

  return(dat)
}

# ==============================================================================
# SUPER LEARNER SETUP
# ==============================================================================

# Configure learners
ranger_learner <- create.Learner("SL.ranger",
                                 params = list(min.node.size = 25),
                                 tune = list(num.trees = c(500, 1000),
                                            mtry = c(3, 4, 5)))

glmnet_learner <- create.Learner("SL.glmnet",
                                 tune = list(alpha = seq(0, 1, .2)))

xgboost_learner <- create.Learner(base_learner = "SL.xgboost",
                                  params = list(minobspernode = 25),
                                  tune = list(max_depth = c(2, 4),
                                             nrounds = c(500, 1000)))

sl_lib <- list("SL.mean", "SL.glm",
               "SL.ranger_1", "SL.ranger_2", "SL.ranger_3",
               "SL.ranger_4", "SL.ranger_5", "SL.ranger_6",
               "SL.xgboost_1", "SL.xgboost_2", "SL.xgboost_3", "SL.xgboost_4",
               "SL.glmnet_1", "SL.glmnet_2", "SL.glmnet_3",
               "SL.glmnet_4", "SL.glmnet_5", "SL.glmnet_6",
               c("SL.glm.interaction", "screen.corP"))

# ==============================================================================
# PARALLELIZATION CONFIGURATION
# ==============================================================================
# Strategy: Parallelize at SCENARIO level, not within CV.SuperLearner
# This prevents nested parallelization which causes memory crashes
#
# Adjust n_parallel_scenarios based on your system:
#   - 2-3 scenarios: Safe for 16-24GB RAM (RECOMMENDED START HERE)
#   - 4 scenarios: Requires 32GB+ RAM
#   - 8 scenarios: Requires 64GB+ RAM
# ==============================================================================

n_parallel_scenarios <- 8  # ADJUST THIS BASED ON YOUR RAM

# Set up future plan for scenario-level parallelization
plan(multisession, workers = n_parallel_scenarios)

# Disable parallelization within CV.SuperLearner to prevent nested parallel
options(mc.cores = 1)

# ==============================================================================
# AIPW SCORING FUNCTION
# ==============================================================================

compute_aipw_scores <- function(dat) {

  # Create exposure and outcome
  exposure <- as.numeric(dat$age_bin)
  outcome <- as.numeric(dat$AFCt)

  # Identify categorical and continuous variables
  categorical_vars <- dat %>%
    select(where(is.factor), -DOR) %>%
    names()

  continuous_vars <- dat %>%
    select(where(is.numeric), -starts_with("imp_"), -age_bin, -AFCt, -age, -DOR) %>%
    names()

  flag_vars <- dat %>% select(starts_with("imp_")) %>% names()

  # Create new dataset
  new_afc <- data.frame(
    outcome,
    exposure,
    dat[, categorical_vars],
    dat[, continuous_vars],
    dat[, flag_vars]
  )

  # Create covariate matrices
  covariates <- c(continuous_vars, categorical_vars, flag_vars)
  covariates_matrix <- data.frame(
    model.matrix(
      formula(paste0("~", paste0(covariates, collapse = "+"))),
      data = new_afc
    )[, -1]
  )
  covariates_matrix_w <- data.frame(covariates_matrix, age_bin = exposure)

  # Set up cross-validation folds
  n <- nrow(new_afc)
  num.folds <- 10
  folds <- sort(seq(n) %% num.folds) + 1
  fold_dat <- tibble(id = 1:n, folds)
  fold_index <- split(fold_dat$id, fold_dat$folds)

  # Fit outcome model (mu)
  cat("    Fitting outcome model (mu)...\n")
  fit_mu <- CV.SuperLearner(
    Y = outcome,
    X = covariates_matrix_w,
    method = "method.NNLS",
    family = gaussian(),
    SL.library = sl_lib,
    cvControl = list(V = num.folds, validRows = fold_index),
    control = list(saveCVFitLibrary = TRUE),
    verbose = FALSE
  )

  # Fit propensity score model (pi)
  cat("    Fitting propensity score model (pi)...\n")
  fit_pi <- CV.SuperLearner(
    Y = exposure,
    X = covariates_matrix,
    method = "method.NNLS",
    family = binomial(),
    SL.library = sl_lib,
    cvControl = list(V = num.folds, validRows = fold_index),
    control = list(saveCVFitLibrary = FALSE),
    verbose = FALSE
  )

  # Get predictions
  pscore <- as.matrix(fit_pi$SL.predict)
  mu_hat <- as.matrix(fit_mu$SL.predict)

  # Get counterfactual predictions
  mu_hat1 <- NULL
  for (i in 1:num.folds) {
    mu_hat1 <- rbind(
      mu_hat1,
      predict(
        fit_mu$AllSL[[i]],
        newdata = transform(covariates_matrix_w[fold_index[[i]], ], age_bin = 1),
        onlySL = TRUE
      )$pred
    )
  }

  mu_hat0 <- NULL
  for (i in 1:num.folds) {
    mu_hat0 <- rbind(
      mu_hat0,
      predict(
        fit_mu$AllSL[[i]],
        newdata = transform(covariates_matrix_w[fold_index[[i]], ], age_bin = 0),
        onlySL = TRUE
      )$pred
    )
  }

  # Compute AIPW scores
  aipw_score <- ((2 * exposure - 1) * (outcome - mu_hat)) /
                ((2 * exposure - 1) * pscore + (1 - exposure)) +
                (mu_hat1 - mu_hat0)

  return(as.numeric(aipw_score))
}

# ==============================================================================
# BEST LINEAR PROJECTION FUNCTION
# ==============================================================================

compute_linear_projection <- function(aipw_scores, edc_data) {

  # Create data frame with DR scores and log-transformed EDCs
  lin_proj <- data.frame(
    dr_scores = aipw_scores,
    edc_data[, env_vars]
  ) %>%
    mutate(across(all_of(env_vars), log))

  # Fit conditional model (all EDCs together)
  model_dr_cond <- lm(dr_scores ~ ., data = lin_proj)
  dr_conditional <- tidy(
    coeftest(model_dr_cond, vcov = vcovHC(model_dr_cond, type = "HC3"))
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(Type = "conditional")

  # Fit unconditional models (one EDC at a time)
  dr_unconditional <- map_dfr(env_vars, function(var) {
    formula_str <- paste("dr_scores ~", var)
    model <- lm(as.formula(formula_str), data = lin_proj)
    tidy(coeftest(model, vcov = vcovHC(model, type = "HC3"))) %>%
      filter(term != "(Intercept)") %>%
      mutate(Type = "unconditional")
  })

  # Combine results
  results <- rbind(dr_unconditional, dr_conditional)

  return(results)
}

# ==============================================================================
# MAIN SENSITIVITY ANALYSIS LOOP (PARALLELIZED)
# ==============================================================================

cat("\n=== Starting Sensitivity Analysis ===\n")
cat("Total combinations: 16\n")
cat("  Outliers: 2 options (in, out)\n")
cat("  EDC trimming: 4 options (90th, 92.5th, 95th, 97.5th percentiles)\n")
cat("  AFC trimming: 2 options (yes, no)\n\n")
cat("Running", n_parallel_scenarios, "scenarios in parallel\n")
cat("Each scenario runs CV.SuperLearner sequentially to avoid memory issues\n\n")

# Define scenario specifications
outlier_options <- c(FALSE, TRUE)
edc_percentiles <- c(0.90, 0.925, 0.95, 0.975)
afc_trim_options <- c(FALSE, TRUE)

# Create grid of all combinations
scenario_grid <- expand.grid(
  outliers_removed = outlier_options,
  edc_percentile = edc_percentiles,
  afc_trimmed = afc_trim_options
)

# Add scenario IDs
scenario_grid$scenario_id <- 1:nrow(scenario_grid)

# Function to process a single scenario
process_scenario <- function(i) {

  scenario <- scenario_grid[i, ]

  # Create log message
  msg <- paste0(
    "\n========================================\n",
    "Scenario ", scenario$scenario_id, " of ", nrow(scenario_grid), "\n",
    "----------------------------------------\n",
    "  Outliers removed: ", scenario$outliers_removed, "\n",
    "  EDC percentile: ", scenario$edc_percentile, "\n",
    "  AFC trimmed: ", scenario$afc_trimmed, "\n",
    "----------------------------------------\n"
  )
  cat(msg)

  # Prepare data for this scenario
  cat("  Preparing data...\n")
  dat <- prepare_data(
    data = afc_base,
    remove_outliers = scenario$outliers_removed,
    edc_percentile = scenario$edc_percentile,
    afc_trim = scenario$afc_trimmed
  )

  # Compute AIPW scores
  cat("  Computing AIPW scores...\n")
  aipw_scores <- compute_aipw_scores(dat)

  # Compute linear projection test statistics
  cat("  Computing linear projection...\n")
  lin_proj_results <- compute_linear_projection(aipw_scores, dat)

  # Add scenario identifiers to results
  lin_proj_results <- lin_proj_results %>%
    mutate(
      scenario_id = scenario$scenario_id,
      outliers_removed = scenario$outliers_removed,
      edc_percentile = scenario$edc_percentile,
      afc_trimmed = scenario$afc_trimmed,
      n_obs = nrow(dat)
    )

  cat("  Completed!\n")

  return(lin_proj_results)
}

# Run scenarios in parallel
cat("Starting parallel execution...\n")

# Explicitly specify globals to export to workers
# This ensures all custom learners and necessary objects are available
all_results <- future_lapply(
  1:nrow(scenario_grid),
  process_scenario,
  future.seed = TRUE,
  future.globals = list(
    scenario_grid = scenario_grid,
    afc_base = afc_base,
    prepare_data = prepare_data,
    compute_aipw_scores = compute_aipw_scores,
    compute_linear_projection = compute_linear_projection,
    env_vars = env_vars,
    the_data = the_data,
    impute_flags = impute_flags,
    outlier_rows = outlier_rows,
    sl_lib = sl_lib,
    # Export the custom learner functions created by create.Learner
    SL.ranger_1 = get("SL.ranger_1", envir = .GlobalEnv),
    SL.ranger_2 = get("SL.ranger_2", envir = .GlobalEnv),
    SL.ranger_3 = get("SL.ranger_3", envir = .GlobalEnv),
    SL.ranger_4 = get("SL.ranger_4", envir = .GlobalEnv),
    SL.ranger_5 = get("SL.ranger_5", envir = .GlobalEnv),
    SL.ranger_6 = get("SL.ranger_6", envir = .GlobalEnv),
    SL.xgboost_1 = get("SL.xgboost_1", envir = .GlobalEnv),
    SL.xgboost_2 = get("SL.xgboost_2", envir = .GlobalEnv),
    SL.xgboost_3 = get("SL.xgboost_3", envir = .GlobalEnv),
    SL.xgboost_4 = get("SL.xgboost_4", envir = .GlobalEnv),
    SL.glmnet_1 = get("SL.glmnet_1", envir = .GlobalEnv),
    SL.glmnet_2 = get("SL.glmnet_2", envir = .GlobalEnv),
    SL.glmnet_3 = get("SL.glmnet_3", envir = .GlobalEnv),
    SL.glmnet_4 = get("SL.glmnet_4", envir = .GlobalEnv),
    SL.glmnet_5 = get("SL.glmnet_5", envir = .GlobalEnv),
    SL.glmnet_6 = get("SL.glmnet_6", envir = .GlobalEnv)
  ),
  future.packages = c("SuperLearner", "tidyverse", "broom", "lmtest",
                      "sandwich", "ranger", "xgboost")
)

cat("\n=== Sensitivity Analysis Complete ===\n\n")

# ==============================================================================
# COMBINE AND SAVE RESULTS
# ==============================================================================

# Combine all results
sensitivity_results <- bind_rows(all_results)

# Save raw results
saveRDS(sensitivity_results, here("output", "sensitivity_analysis_results.rds"))

cat("Results saved to: output/sensitivity_analysis_results.rds\n")

# ==============================================================================
# CREATE VISUALIZATION: test_scatter.png
# ==============================================================================

cat("\nCreating test_scatter.png visualization...\n")

# Create descriptive labels for scenarios
plot_data <- sensitivity_results %>%
  mutate(
    outlier_label = ifelse(outliers_removed, "Out", "In"),
    edc_label = paste0(edc_percentile * 100, "th"),
    afc_label = ifelse(afc_trimmed, "AFC≤30", "AFC orig"),
    scenario_label = paste0(
      "Outliers: ", outlier_label,
      " | EDC: ", edc_label,
      " | ", afc_label
    )
  )


save(plot_data, 
     file = here("data", "sensitivity_data.Rdata"))

table_data <- readRDS(here("output", "linear_projection_output.rds"))

# Create plot of DR learner test statistics
ggplot() +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(data = plot_data,
             aes(x = statistic, y = reorder(term, statistic), shape = Type),
             size = 3, color = "gray", alpha = 0.5) +
  geom_point(data = table_data,
             aes(x = statistic, y = reorder(term, statistic), shape = Type),
             size = 3, color = "black") +
  scale_shape_manual(values = c(16, 3)) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0))

ggsave(filename = here("figures", "teststat_scatter_sensitivity.png"),
       width = 16, height = 16, units = "cm", dpi = 300)

cat("\n=== Analysis Complete ===\n")

# ==============================================================================
# CLEANUP AND MEMORY MANAGEMENT
# ==============================================================================

# Close parallel workers to free memory
plan(sequential)

cat("\nMemory cleanup complete. Parallel workers closed.\n")
