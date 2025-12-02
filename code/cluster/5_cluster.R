#!/usr/bin/env Rscript
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)
# ==============================================================================
# CLUSTER-OPTIMIZED SENSITIVITY ANALYSIS FOR AFC PROJECT
# ==============================================================================
# This script is designed to run on an HPC cluster with SLURM.
# It processes a subset of 24 sensitivity analysis scenarios based on
# the SLURM_ARRAY_TASK_ID passed as a command line argument.
#
# Called by: earth.sh (SLURM submission script)
# ==============================================================================

# Get array task ID from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No array task ID provided. This script should be called with a task ID argument.")
}
task_id <- as.integer(args[1])

cat("\n")
cat("================================================================================\n")
cat("AFC SENSITIVITY ANALYSIS - CLUSTER EXECUTION\n")
cat("================================================================================\n")
cat("Task ID:", task_id, "\n")
cat("Started at:", as.character(Sys.time()), "\n")
cat("Node:", Sys.info()["nodename"], "\n")
cat("================================================================================\n\n")

# ==============================================================================
# LOAD PACKAGES
# ==============================================================================

cat("Loading required packages...\n")

# Install missForest dependencies explicitly with better error reporting
cat("Installing missForest dependencies...\n")

# Install foreach (required by doRNG)
if (!requireNamespace("foreach", quietly = TRUE)) {
  cat("  Installing foreach...\n")
  install.packages("foreach", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

# Install rngtools (required by doRNG)
if (!requireNamespace("rngtools", quietly = TRUE)) {
  cat("  Installing rngtools...\n")
  install.packages("rngtools", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

# Install iterators (required by doRNG)
if (!requireNamespace("iterators", quietly = TRUE)) {
  cat("  Installing iterators...\n")
  install.packages("iterators", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

# Install doRNG (required by missForest)
if (!requireNamespace("doRNG", quietly = TRUE)) {
  cat("  Installing doRNG...\n")
  install.packages("doRNG", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

# Install randomForest (required by missForest)
if (!requireNamespace("randomForest", quietly = TRUE)) {
  cat("  Installing randomForest...\n")
  install.packages("randomForest", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

# Install missForest
if (!requireNamespace("missForest", quietly = TRUE)) {
  cat("  Installing missForest...\n")
  install.packages("missForest", repos = "https://cloud.r-project.org/", lib = userLib,
                   dependencies = TRUE, type = "source")
}

cat("Dependency installation complete.\n\n")

# Verify library paths
cat("Current library paths:\n")
print(.libPaths())
cat("\n")

# Now load all packages with pacman
pacman::p_load(
  rio,
  here,
  tidyverse,
  parallel,
  SuperLearner,
  broom,
  lmtest,
  sandwich,
  ranger,
  xgboost,
  missForest,
  doParallel
)
cat("Packages loaded successfully.\n")

# Verify missForest is loaded - TERMINATE if not
if ("missForest" %in% loadedNamespaces()) {
  cat("✓ missForest loaded successfully\n\n")
} else {
  cat("\n")
  cat("================================================================================\n")
  cat("ERROR: missForest failed to load properly\n")
  cat("================================================================================\n")
  cat("This script requires missForest for imputation.\n")
  cat("Please check the error messages above to diagnose the installation failure.\n")
  cat("Common issues:\n")
  cat("  - Missing Fortran compiler (needed for randomForest)\n")
  cat("  - Missing system dependencies\n")
  cat("  - Insufficient permissions in user library\n")
  cat("================================================================================\n")
  cat("TERMINATING EARLY to save compute time.\n")
  cat("================================================================================\n")
  quit(status = 1, save = "no")
}
cat("\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

set.seed(123)

# Determine which scenarios this task should process
# With 6 jobs and 24 scenarios: each job processes 4 scenarios
scenarios_per_job <- 4
scenario_start <- (task_id - 1) * scenarios_per_job + 1
scenario_end <- task_id * scenarios_per_job

cat("This task will process scenarios", scenario_start, "to", scenario_end, "\n\n")

# ==============================================================================
# LOAD DATA AND DEFINE VARIABLES
# ==============================================================================

cat("Loading base data...\n")
load(here("data", "imputed_EARTH.Rdata"))
afc_base <- imputed_data %>% rename(AFCt = AFC_t)

# Define variables
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "BPA", "BP", "MP", "PP", "Hg")

the_data <- c("AFCt", "year", "month", "DOR", "eversmok",
              "age", "bmi", "races", "educ1", "smokstat", "previousIVF",
              "previousIUI", "gravid", env_vars)

impute_flags <- c("imp_sgratio_pht", "imp_smokstat", "imp_races", "imp_bmi",
                  "imp_previousIUI", "imp_age", "imp_educ1", "imp_Hg",
                  "imp_mBP", "imp_mCNP", "imp_B_PB", "imp_AFScanDate")

cat("Data loaded successfully.\n")
cat("Sample size:", nrow(afc_base), "\n")
cat("Number of EDCs:", length(env_vars), "\n\n")

# ==============================================================================
# DATA PREPARATION FUNCTION
# ==============================================================================

prepare_data <- function(data,
                        age_threshold = 35,
                        afc_impute = FALSE,
                        edc_impute = FALSE) {

  # Start with base data
  dat <- data[, c(the_data, impute_flags)]

  # Set AFC>30 to missing and re-impute if specified
  if (afc_impute) {
    n_afc_missing <- sum(dat$AFCt > 30)
    dat$AFCt[dat$AFCt > 30] <- NA
    cat("  Set AFC>30 to missing (", n_afc_missing, "obs)\n")
  }

  # Set EDCs>97.5th percentile to missing and re-impute if specified
  if (edc_impute) {
    n_missing_total <- 0
    for (var in env_vars) {
      threshold <- quantile(dat[[var]], 0.975)
      n_missing <- sum(dat[[var]] > threshold)
      dat[[var]][dat[[var]] > threshold] <- NA
      n_missing_total <- n_missing_total + n_missing
    }
    cat("  Set", n_missing_total, "EDC values >97.5th percentile to missing\n")
  }

  # Re-impute using missForest if either AFC or EDC imputation is TRUE
  if (afc_impute || edc_impute) {
    cat("  Running missForest imputation...\n")
    # Prepare data for imputation (exclude imputation flags)
    data_for_imputation <- dat %>% select(-starts_with("imp_"))

    # Run missForest
    imputation_result <- missForest(
      xmis = data.frame(data_for_imputation),
      maxiter = 50,
      ntree = 2000,
      verbose = FALSE,
      parallelize = "no"
    )

    # Extract imputed data and merge back with imputation flags
    dat_imputed <- imputation_result$ximp
    dat <- bind_cols(dat_imputed, dat %>% select(starts_with("imp_")))

    cat("  missForest imputation complete (OOB error:",
        round(imputation_result$OOBerror, 4), ")\n")
  }

  # Create age binary variable with specified threshold
  dat$age_bin <- as.numeric(dat$age >= age_threshold)
  cat("  Using age threshold:", age_threshold, "years\n")

  return(dat)
}

# ==============================================================================
# SUPER LEARNER SETUP
# ==============================================================================

cat("Setting up SuperLearner libraries...\n")

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

cat("SuperLearner libraries configured.\n\n")

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
# CREATE SCENARIO GRID
# ==============================================================================

# Define scenario specifications
age_thresholds <- c(35, 36, 37, 38, 39, 40)
afc_impute_options <- c(FALSE, TRUE)
edc_impute_options <- c(FALSE, TRUE)

# Create grid of all combinations
scenario_grid <- expand.grid(
  age_threshold = age_thresholds,
  afc_impute = afc_impute_options,
  edc_impute = edc_impute_options
)

# Add scenario IDs
scenario_grid$scenario_id <- 1:nrow(scenario_grid)

cat("Total scenarios in grid:", nrow(scenario_grid), "\n\n")

# ==============================================================================
# PROCESS SCENARIOS FOR THIS TASK
# ==============================================================================

# Filter to scenarios for this task
scenarios_to_run <- scenario_grid %>%
  filter(scenario_id >= scenario_start, scenario_id <= scenario_end)

cat("================================================================================\n")
cat("PROCESSING", nrow(scenarios_to_run), "SCENARIOS\n")
cat("================================================================================\n\n")

# Storage for results
all_results <- list()

# Process each scenario sequentially
for (i in 1:nrow(scenarios_to_run)) {

  scenario <- scenarios_to_run[i, ]

  # Create log message
  msg <- paste0(
    "\n========================================\n",
    "Scenario ", scenario$scenario_id, " of 24\n",
    "----------------------------------------\n",
    "  Age threshold: ", scenario$age_threshold, " years\n",
    "  AFC>30 imputation: ", scenario$afc_impute, "\n",
    "  EDC>97.5th imputation: ", scenario$edc_impute, "\n",
    "----------------------------------------\n"
  )
  cat(msg)

  # Prepare data for this scenario
  cat("  Preparing data...\n")
  dat <- prepare_data(
    data = afc_base,
    age_threshold = scenario$age_threshold,
    afc_impute = scenario$afc_impute,
    edc_impute = scenario$edc_impute
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
      age_threshold = scenario$age_threshold,
      afc_impute = scenario$afc_impute,
      edc_impute = scenario$edc_impute,
      n_obs = nrow(dat)
    )

  # Store results
  all_results[[i]] <- lin_proj_results

  cat("  Completed!\n")
}

# ==============================================================================
# COMBINE AND SAVE RESULTS
# ==============================================================================

cat("\n================================================================================\n")
cat("COMBINING AND SAVING RESULTS\n")
cat("================================================================================\n\n")

# Combine all results for this task
task_results <- bind_rows(all_results)

# Save results for this task
output_file <- here("output", paste0("sensitivity_task_", task_id, ".rds"))
saveRDS(task_results, file = output_file)

cat("Results saved to:", output_file, "\n")
cat("  Rows:", nrow(task_results), "\n")
cat("  Scenarios:", n_distinct(task_results$scenario_id), "\n")

cat("\n================================================================================\n")
cat("TASK COMPLETED SUCCESSFULLY\n")
cat("Finished at:", as.character(Sys.time()), "\n")
cat("================================================================================\n\n")
