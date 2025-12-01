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
  xgboost
)

# ==============================================================================
# SENSITIVITY ANALYSIS: Outliers Removed
# ==============================================================================
# This program re-runs the main analysis with consensus outliers removed
# (no trimming or resetting of extreme values)
#
# Output: Best linear projection test statistics for plotting as blue dots
# ==============================================================================

set.seed(123)

cat("\n=== Sensitivity Analysis: Outliers Removed ===\n")
cat("This analysis removes consensus outliers and re-runs the full pipeline\n")
cat("No trimming or resetting of extreme AFC or EDC values\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load base data
load(here("data", "imputed_EARTH.Rdata"))
afc_base <- imputed_data %>% rename(AFCt = AFC_t)

# Load outlier IDs (flagged by 2+ methods from program 2)
outliers_id <- readRDS(here("output", "outliers_id.rds"))
outlier_rows <- outliers_id$row_id

cat("Total observations:", nrow(afc_base), "\n")
cat("Outliers to remove:", length(outlier_rows), "\n")

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
# PREPARE DATA (REMOVE OUTLIERS ONLY)
# ==============================================================================

cat("\nPreparing data...\n")

# Start with base data
dat <- afc_base[, c(the_data, impute_flags)]

# Store row IDs before filtering
dat$row_id <- as.numeric(rownames(dat))

# Remove outliers
dat <- dat %>% filter(!row_id %in% outlier_rows)

cat("Observations after removing outliers:", nrow(dat), "\n")

# Create age binary variable
dat$age_bin <- as.numeric(dat$age >= 35)

# Remove row_id before analysis
dat <- dat %>% select(-row_id)

# ==============================================================================
# SUPER LEARNER SETUP
# ==============================================================================

cat("\nConfiguring SuperLearner...\n")

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

# Set parallel processing within CV.SuperLearner for efficiency
options(mc.cores = parallel::detectCores() - 2)

cat("Using", getOption("mc.cores"), "cores for parallel processing\n")

# ==============================================================================
# COMPUTE AIPW SCORES
# ==============================================================================

cat("\n=== Computing AIPW Scores ===\n")

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
cat("Fitting outcome model (mu)...\n")
cat("  This may take several minutes...\n")

fit_mu <- CV.SuperLearner(
  Y = outcome,
  X = covariates_matrix_w,
  method = "method.NNLS",
  family = gaussian(),
  SL.library = sl_lib,
  cvControl = list(V = num.folds, validRows = fold_index),
  control = list(saveCVFitLibrary = TRUE),
  parallel = "multicore",
  verbose = FALSE
)

cat("Outcome model complete.\n")

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
  parallel = "multicore",
  verbose = FALSE
)

cat("Propensity score model complete.\n")

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

# ==============================================================================
# BEST LINEAR PROJECTION
# ==============================================================================

cat("\n=== Computing Best Linear Projection ===\n")

# Create data frame with DR scores and log-transformed EDCs
lin_proj <- data.frame(
  dr_scores = as.numeric(aipw_score),
  dat[, env_vars]
) %>%
  mutate(across(all_of(env_vars), log))

# Fit conditional model (all EDCs together)
cat("Fitting conditional model...\n")
model_dr_cond <- lm(dr_scores ~ ., data = lin_proj)
dr_conditional <- tidy(
  coeftest(model_dr_cond, vcov = vcovHC(model_dr_cond, type = "HC3"))
) %>%
  filter(term != "(Intercept)") %>%
  mutate(Type = "conditional")

# Fit unconditional models (one EDC at a time)
cat("Fitting unconditional models...\n")
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
    analysis = "outliers_removed",
    n_obs = nrow(dat)
  )

cat("Best linear projection complete.\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n=== Saving Results ===\n")

# Save detailed results
saveRDS(results, here("output", "sensitivity_outliers_removed.rds"))
cat("Detailed results saved to: output/sensitivity_outliers_removed.rds\n")

# Save AIPW scores for potential further analysis
saveRDS(aipw_score, here("output", "aipw_scores_outliers_removed.rds"))
cat("AIPW scores saved to: output/aipw_scores_outliers_removed.rds\n")

# Display summary
cat("\n=== Summary of Results ===\n")
cat("Sample size:", nrow(dat), "\n")
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
cat("These results can be added as blue dots to the test statistic scatter plot\n")
