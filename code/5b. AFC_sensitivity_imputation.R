pacman::p_load(
  rio,
  here,
  tidyverse,
  parallel,
  doParallel,
  SuperLearner,
  broom,
  lmtest,
  sandwich,
  ranger,
  xgboost,
  missForest,
  dbscan,
  MVN
)

# ==============================================================================
# SENSITIVITY ANALYSIS: Extreme Values Imputed
# ==============================================================================
# This program sets extreme AFC and EDC values to missing, re-imputes them,
# re-identifies outliers, and re-runs the full analysis
#
# Steps:
# 1. Set AFC > 30 and EDC > 92.5th percentile to NA
# 2. Impute using missForest
# 3. Re-run outlier detection (4 methods: Mahalanobis, PCA, LOF, Cook's D)
# 4. Run SuperLearner + Best Linear Projection (WITHOUT removing new outliers)
# 5. Bootstrap option for uncertainty quantification
#
# Output: Best linear projection test statistics for plotting as gray dots
# ==============================================================================

set.seed(123)

cat("\n=== Sensitivity Analysis: Extreme Values Imputed ===\n")
cat("Sets AFC>30 and EDC>92.5th percentile to missing\n")
cat("Re-imputes, re-detects outliers, and re-runs analysis\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set to TRUE to run bootstrap, FALSE for single analysis
RUN_BOOTSTRAP <- T
N_BOOTSTRAP <- 50  # Number of bootstrap iterations (if RUN_BOOTSTRAP = TRUE)

cat("Bootstrap mode:", ifelse(RUN_BOOTSTRAP, "ENABLED", "DISABLED"), "\n")
if (RUN_BOOTSTRAP) {
  cat("Bootstrap iterations:", N_BOOTSTRAP, "\n")
}
cat("\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load base data
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

cat("Total observations:", nrow(afc_base), "\n")

# ==============================================================================
# STEP 1: SET EXTREME VALUES TO MISSING
# ==============================================================================

cat("\n=== Step 1: Setting Extreme Values to Missing ===\n")

dat <- afc_base[, c(the_data, impute_flags)]

# Store original values for comparison
n_missing_before <- sum(is.na(dat))

# Set AFC > 30 to NA
afc_threshold <- 30
n_afc_extreme <- sum(dat$AFCt > afc_threshold, na.rm = TRUE)
dat$AFCt[dat$AFCt > afc_threshold] <- NA

cat("AFC > 30 set to missing:", n_afc_extreme, "\n")

# Set EDC values > 97.5th percentile to NA
edc_percentile <- 0.975
n_edc_extreme <- 0

for (var in env_vars) {
  threshold <- quantile(dat[[var]], edc_percentile, na.rm = TRUE)
  n_extreme <- sum(dat[[var]] > threshold, na.rm = TRUE)
  dat[[var]][dat[[var]] > threshold] <- NA
  n_edc_extreme <- n_edc_extreme + n_extreme
  cat("  ", var, "> ", round(threshold, 3), " set to missing: ", n_extreme, "\n", sep = "")
}

n_missing_after <- sum(is.na(dat))
cat("\nTotal missing values before:", n_missing_before, "\n")
cat("Total missing values after:", n_missing_after, "\n")
cat("New missing values created:", n_missing_after - n_missing_before, "\n")

# ==============================================================================
# STEP 2: RE-IMPUTE USING MISSFOREST
# ==============================================================================

cat("\n=== Step 2: Re-imputing with missForest ===\n")

# Prepare data for imputation (only numeric and factor variables)
impute_data <- dat %>%
  select(AFCt, all_of(env_vars), age, bmi, year, month,
         races, educ1, smokstat, previousIVF, previousIUI, gravid)

cat("Starting missForest imputation...\n")
cat("  This may take 5-15 minutes...\n")

n_cores <- detectCores() - 2
registerDoParallel(cores = n_cores)
# Run missForest with optimized settings for speed
impute_result <- missForest(
  xmis = data.frame(impute_data),
  maxiter = 50,        # Limit iterations for speed
  ntree = 2000,        # Fewer trees for speed
  parallelize = "forests",  # Parallelize across trees
  verbose = TRUE
)
stopImplicitCluster()

cat("\nImputation complete!\n")
cat("OOB error (NRMSE):", round(impute_result$OOBerror[1], 4), "\n")
cat("OOB error (PFC):", round(impute_result$OOBerror[2], 4), "\n")

# Replace imputed values back into dataset
dat_imputed <- dat
dat_imputed[, names(impute_data)] <- impute_result$ximp

# Save imputed dataset
save(dat_imputed, file = here("data", "imputed_EARTH_5b.Rdata"))
cat("Imputed data saved to: data/imputed_EARTH_5b.Rdata\n")

# ==============================================================================
# STEP 3: RE-IDENTIFY OUTLIERS
# ==============================================================================

cat("\n=== Step 3: Re-identifying Outliers ===\n")

# Prepare data for multivariate outlier detection
multivar_data <- dat_imputed %>%
  select(AFCt, age, all_of(env_vars)) %>%
  mutate(row_id = row_number())

cat("Sample size for outlier detection:", nrow(multivar_data), "\n")

# METHOD 1: Mahalanobis Distance
cat("\n--- Method 1: Mahalanobis Distance ---\n")

mahala_vars <- multivar_data %>% select(-row_id)
mahala_center <- colMeans(mahala_vars)
mahala_cov <- cov(mahala_vars)
mahala_dist <- mahalanobis(mahala_vars, mahala_center, mahala_cov)

# Chi-square critical value at p < 0.001
df <- ncol(mahala_vars)
crit_value <- qchisq(0.999, df)

multivar_data$mahala_dist <- mahala_dist
multivar_data$mahala_outlier <- mahala_dist > crit_value

cat("Mahalanobis outliers (p < 0.001):", sum(multivar_data$mahala_outlier), "\n")

# METHOD 2: PCA-Based Outliers
cat("\n--- Method 2: PCA-Based Outliers ---\n")

pca_result <- prcomp(mahala_vars, scale. = TRUE, center = TRUE)
pca_scores <- as.data.frame(pca_result$x)
n_pcs <- min(5, ncol(pca_scores))
pca_scores_subset <- pca_scores[, 1:n_pcs]

pc_dist <- sqrt(rowSums(pca_scores_subset^2))
pc_threshold <- quantile(pc_dist, 0.99)

multivar_data$pc_dist <- pc_dist
multivar_data$pc_outlier <- pc_dist > pc_threshold

cat("PCA outliers (top 1%):", sum(multivar_data$pc_outlier), "\n")

# METHOD 3: Local Outlier Factor (LOF)
cat("\n--- Method 3: Local Outlier Factor ---\n")

lof_scores <- lof(mahala_vars, k = 10)
lof_threshold <- 1.5

multivar_data$lof_score <- lof_scores
multivar_data$lof_outlier <- lof_scores > lof_threshold

cat("LOF outliers (LOF > 1.5):", sum(multivar_data$lof_outlier), "\n")

# METHOD 4: Cook's Distance
cat("\n--- Method 4: Cook's Distance ---\n")

cooksd_model <- lm(AFCt ~ ., data = mahala_vars)
cooksd <- cooks.distance(cooksd_model)
cooksd_threshold <- 4 / nrow(mahala_vars)

multivar_data$cooksd <- cooksd
multivar_data$high_cooksd <- cooksd > cooksd_threshold

cat("High Cook's D observations:", sum(multivar_data$high_cooksd), "\n")

# Consensus: Count how many methods flag each observation
multivar_data$n_flags <- rowSums(multivar_data[, c("mahala_outlier",
                                                     "pc_outlier",
                                                     "lof_outlier",
                                                     "high_cooksd")])

consensus_outliers <- multivar_data %>% filter(n_flags >= 2)

cat("\n=== OUTLIER SUMMARY ===\n")
cat("Observations flagged by 2+ methods:", nrow(consensus_outliers), "\n")

# Save outlier results
outliers_5b <- consensus_outliers %>%
  select(row_id, AFCt, age, n_flags,
         mahala_outlier, pc_outlier, lof_outlier, high_cooksd)

saveRDS(outliers_5b, here("output", "outliers_5b.rds"))
write_csv(multivar_data, here("misc", "outlier_analysis_5b.csv"))

cat("Outlier results saved to: output/outliers_5b.rds\n")
cat("Full outlier analysis saved to: misc/outlier_analysis_5b.csv\n")

# ==============================================================================
# STEP 4: PREPARE DATA FOR ANALYSIS
# ==============================================================================

cat("\n=== Step 4: Preparing Data for Analysis ===\n")

# Use imputed data (do NOT remove newly identified outliers)
dat_final <- dat_imputed[, c(the_data, impute_flags)]

# Create age binary variable
dat_final$age_bin <- as.numeric(dat_final$age >= 35)

cat("Final sample size:", nrow(dat_final), "\n")
cat("Note: Newly identified outliers are NOT removed for this analysis\n")

# ==============================================================================
# STEP 5: SUPER LEARNER SETUP
# ==============================================================================

cat("\n=== Step 5: Configuring SuperLearner ===\n")

# Configure learners (same as main analysis)
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

# Set parallel processing
options(mc.cores = parallel::detectCores() - 2)
cat("Using", getOption("mc.cores"), "cores for parallel processing\n")

# ==============================================================================
# FUNCTION: COMPUTE ANALYSIS
# ==============================================================================

compute_analysis <- function(data, iteration = NULL) {

  if (!is.null(iteration)) {
    cat("\n--- Bootstrap Iteration", iteration, "---\n")
  } else {
    cat("\n=== Computing AIPW Scores ===\n")
  }

  # Create exposure and outcome
  exposure <- as.numeric(data$age_bin)
  outcome <- as.numeric(data$AFCt)

  # Identify categorical and continuous variables
  categorical_vars <- data %>%
    select(where(is.factor), -DOR) %>%
    names()

  continuous_vars <- data %>%
    select(where(is.numeric), -starts_with("imp_"), -age_bin, -AFCt, -age, -DOR) %>%
    names()

  flag_vars <- data %>% select(starts_with("imp_")) %>% names()

  # Create new dataset
  new_afc <- data.frame(
    outcome,
    exposure,
    data[, categorical_vars],
    data[, continuous_vars],
    data[, flag_vars]
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
  cat("Fitting outcome model (mu)...\n")

  fit_mu <- CV.SuperLearner(
    Y = outcome,
    X = covariates_matrix_w,
    method = "method.NNLS",
    family = gaussian(),
    SL.library = sl_lib,
    cvControl = list(V = num.folds, validRows = fold_index),
    control = list(saveCVFitLibrary = TRUE),
    parallel = "seq",
    verbose = FALSE
  )

  # Fit propensity score model (pi)
  cat("Fitting propensity score model (pi)...\n")

  fit_pi <- CV.SuperLearner(
    Y = exposure,
    X = covariates_matrix,
    method = "method.NNLS",
    family = binomial(),
    SL.library = sl_lib,
    cvControl = list(V = num.folds, validRows = fold_index),
    control = list(saveCVFitLibrary = FALSE),
    parallel = "seq",
    verbose = FALSE
  )

  # Get predictions
  pscore <- as.matrix(fit_pi$SL.predict)
  mu_hat <- as.matrix(fit_mu$SL.predict)

  # Get counterfactual predictions
  cat("Computing counterfactual predictions...\n")

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

  cat("AIPW scores computed.\n")

  # Best Linear Projection
  cat("Computing best linear projection...\n")

  # Create data frame with DR scores and log-transformed EDCs
  lin_proj <- data.frame(
    dr_scores = as.numeric(aipw_score),
    data[, env_vars]
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
  results <- rbind(dr_unconditional, dr_conditional) %>%
    mutate(
      analysis = "imputation",
      n_obs = nrow(data),
      n_outliers_identified = nrow(consensus_outliers)
    )

  if (!is.null(iteration)) {
    results <- results %>% mutate(bootstrap_iter = iteration)
  }

  cat("Analysis complete.\n")

  return(list(results = results, aipw_scores = aipw_score))
}

# ==============================================================================
# STEP 6: RUN ANALYSIS
# ==============================================================================

if (!RUN_BOOTSTRAP) {

  # Single analysis without bootstrap
  cat("\n=== Running Single Analysis ===\n")

  analysis_output <- compute_analysis(dat_final)
  results <- analysis_output$results
  aipw_scores <- analysis_output$aipw_scores

  # Save results
  saveRDS(results, here("output", "sensitivity_imputation.rds"))
  saveRDS(aipw_scores, here("output", "aipw_scores_imputation.rds"))

  cat("\nResults saved to: output/sensitivity_imputation.rds\n")
  cat("AIPW scores saved to: output/aipw_scores_imputation.rds\n")

} else {

  # Bootstrap analysis
  cat("\n=== Running Bootstrap Analysis ===\n")
  cat("This will take a long time...\n\n")

  # First run main analysis
  cat("Running main analysis...\n")
  main_output <- compute_analysis(dat_final)
  main_results <- main_output$results

  # Initialize bootstrap results storage
  bootstrap_results <- list()
  bootstrap_results[[1]] <- main_results %>% mutate(bootstrap_iter = 0)

  # Run bootstrap iterations
  for (b in 1:N_BOOTSTRAP) {

    # Bootstrap sample
    boot_indices <- sample(1:nrow(dat_final), replace = TRUE)
    dat_boot <- dat_final[boot_indices, ]

    # Run analysis on bootstrap sample
    boot_output <- compute_analysis(dat_boot, iteration = b)
    bootstrap_results[[b + 1]] <- boot_output$results
  }

  # Combine all bootstrap results
  all_bootstrap_results <- bind_rows(bootstrap_results)

  # Save bootstrap results (individual test statistics for each iteration)
  saveRDS(all_bootstrap_results, here("output", "sensitivity_imputation_bootstrap.rds"))

  cat("\nBootstrap results saved to: output/sensitivity_imputation_bootstrap.rds\n")

  # Use main results for plotting
  results <- main_results
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== Summary of Results ===\n")
cat("Sample size:", unique(results$n_obs), "\n")
cat("Outliers re-identified:", unique(results$n_outliers_identified), "\n")
cat("Number of EDCs:", length(env_vars), "\n")
cat("Total test statistics:", nrow(results), "\n")

summary_stats <- results %>%
  group_by(Type) %>%
  summarise(
    n_sig_10 = sum(abs(statistic) > 1.645),
    n_sig_05 = sum(abs(statistic) > 1.96),
    n_sig_01 = sum(abs(statistic) > 2.576),
    mean_abs_stat = mean(abs(statistic)),
    max_abs_stat = max(abs(statistic)),
    .groups = "drop"
  )

print(summary_stats)

cat("\n=== Analysis Complete ===\n")
cat("These results can be added as gray dots to the test statistic scatter plot\n")
