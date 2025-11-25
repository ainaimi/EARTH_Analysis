pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  parallel,
  scales,
  haven,
  VIM,
  naniar,
  missForest,
  doParallel,
  foreach
)

# set the theme for figures
thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    text = element_text(size = 20)
  )
theme_set(thm)

a <- read_sas(here("data", "afc_mixtures.sas7bdat"))

## step 1: select relevant variables to construct SG versions of EDCs

a_ <- a %>% select(
  id,              # Keep original ID for crosswalk
  sgratio_pht,
  sgratio_bpa,
  sgratio_opfr,

  AFC_t,
  smokstat,
  races,
  bmi,
  previousIUI,
  previousIVF,
  age,
  eversmok,
  educ1,
  gravid,
  Hg,
  mBP,
  miBP,
  mCNP,
  mCOP,
  mECPP,
  mEHHP,
  mEHP,
  mEOHP,
  mCPP,
  mEP1,
  mBZP1,
  dehp,
  BPA,
  BPS,
  BPF,
  BP_3,
  B_PB,
  M_PB,
  P_PB,
  TCS,
  BDCIPP,
  DPHP,
  ipPPP,
  sartnew2,
  AFScanDate
)

skimr::skim(a_)

#### Missing Data Analysis ####

# 1. Examine missing data patterns
# Count missing values per variable
missing_summary <- a_ %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = (n_missing / nrow(a_)) * 100,
    n_complete = nrow(a_) - n_missing
  ) %>%
  arrange(desc(pct_missing))

print(missing_summary, n = 40)

# Save missing data summary
saveRDS(missing_summary, file = here("output", "missing_data_summary.rds"))

# 2. Visualize missing data patterns
# Aggregation plot showing combinations of missingness
aggr_plot <- VIM::aggr(a_,
                       col = c('navyblue', 'red'),
                       numbers = TRUE,
                       sortVars = TRUE,
                       labels = names(a_),
                       cex.axis = 0.7,
                       gap = 3,
                       ylab = c("Histogram of missing data", "Pattern"))

# 3. Remove variables with > 40% missing
remov_var <- missing_summary %>% 
  filter(pct_missing > 40) %>% 
  select(variable)

#### because sgratio_opfr is missing >40, we also need to remove 
#### BDCIPP, DPHP, and ipPPP
a_ <- a_ %>% select(-remov_var$variable, -BDCIPP, -DPHP, -ipPPP)

# Additional missingness analysis
# Check if missingness is related to observed variables
miss_var_summary <- miss_var_summary(a_)
print(miss_var_summary, n = 40)

# Visualize missingness by variable
miss_var_plot <- ggplot(miss_var_summary,
                        aes(x = reorder(variable, pct_miss), y = pct_miss)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Variable",
       y = "Percent of Missing Values",
       title = "Missing Data by Variable") +
  theme_classic(base_size = 12)

ggsave(here("figures", "missing_by_variable.png"),
       miss_var_plot, width = 10, height = 8, dpi = 300)


#### Random Forest Imputation ####

# Prepare data for imputation
# Convert character variables to factors if needed
a_for_imputation <- a_ %>%
  mutate(across(c(smokstat, races, educ1, sartnew2, 
                  previousIUI, previousIVF, eversmok, 
                  gravid), as.factor),
         year = year(AFScanDate),
         month = month(AFScanDate),
         sartnew2 = if_else(sartnew2=="", NA, sartnew2)) %>% 
  select(-AFScanDate, -id)

# Perform random forest imputation
# ntree: number of trees in random forest
# mtry: number of variables randomly sampled at each split
# verbose: show progress
n_cores <- detectCores() - 2
registerDoParallel(cores = n_cores)
set.seed(123)  # For reproducibility
imputation_result <- missForest(
  xmis = data.frame(a_for_imputation),
  maxiter = 50,
  ntree = 2000,
  verbose = TRUE,
  parallelize = "variables"
)
stopImplicitCluster()

# Extract imputed data
a_imputed <- imputation_result$ximp

# Extract imputation error estimates
# OOB error for continuous variables (NRMSE: normalized root mean squared error)
# OOB error for categorical variables (proportion of falsely classified)
imputation_error <- data.frame(
  NRMSE = imputation_result$OOBerror[1],  # For continuous variables
  PFC = imputation_result$OOBerror[2]      # For categorical variables
)

print("Imputation Out-of-Bag Errors:")
print(imputation_error)

# Save imputation errors
saveRDS(imputation_error, file = here("output", "imputation_error.rds"))

# 5. Compare distributions before and after imputation
# Select 12 variables with the most missingness
top_missing_vars <- missing_summary %>%
  filter(n_missing > 0, variable %in% names(a_for_imputation)) %>%
  arrange(desc(n_missing)) %>%
  dplyr::slice(1:12) %>%
  pull(variable)

# Prepare data for ggplot
comparison_data <- map_dfr(top_missing_vars, function(var) {
  if (is.numeric(a_for_imputation[[var]])) {
    # Create density data for original (observed values only)
    dens_orig <- density(a_for_imputation[[var]], na.rm = TRUE)
    orig_df <- data.frame(
      x = dens_orig$x,
      y = dens_orig$y,
      variable = var,
      type = "Original (observed)"
    )

    # Create density data for imputed (all values)
    dens_imp <- density(a_imputed[[var]], na.rm = TRUE)
    imp_df <- data.frame(
      x = dens_imp$x,
      y = dens_imp$y,
      variable = var,
      type = "Imputed (complete)"
    )

    bind_rows(orig_df, imp_df)
  } else {
    NULL
  }
})

# Create ggplot
imputation_comparison_plot <- ggplot(comparison_data,
                                     aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  scale_color_manual(values = c("Original (observed)" = "blue",
                                 "Imputed (complete)" = "red")) +
  scale_linetype_manual(values = c("Original (observed)" = "solid",
                                    "Imputed (complete)" = "solid")) +
  labs(x = "Value",
       y = "Density",
       color = NULL,
       linetype = NULL,
       title = "Distribution Comparison: Original vs Imputed Data",
       subtitle = "12 variables with most missing data") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "gray90"))

imputation_comparison_plot

ggsave(here("figures", "imputation_comparison.png"),
       imputation_comparison_plot, width = 12, height = 10, dpi = 300)

# 6. Create unique identifiers and crosswalk
# Generate new sequential study ID
study_id <- 1:nrow(a_)

# Create crosswalk linking original ID to new unique ID
id_crosswalk <- data.frame(
  study_id = study_id,
  original_id = a_$id
)

# Save crosswalk securely
saveRDS(id_crosswalk, file = here("data", "id_crosswalk.rds"))

# 6b. Create indicator variables for imputed values
# This allows checking if imputation affects results in sensitivity analyses
imputation_indicators <- a_ %>%
  select(-id) %>%  # Remove original ID
  mutate(across(everything(), ~as.integer(is.na(.)), .names = "imp_{.col}")) %>%
  mutate(study_id = study_id, .before = 1) %>%
  select(starts_with("imp_"), study_id) %>%
  # Remove indicators that are all zeros (no missingness)
  select(where(~!all(. == 0)), study_id) %>%
  # Remove linearly dependent imputation flags (keep one representative from each group)
  select(-imp_sgratio_bpa,        # identical to imp_sgratio_pht
         -imp_miBP, -imp_mECPP, -imp_mEHHP, -imp_mEHP, -imp_mEOHP,  # identical to imp_sgratio_pht
         -imp_mCPP, -imp_mEP1, -imp_mBZP1, -imp_dehp, -imp_BPA,     # identical to imp_sgratio_pht
         -imp_mCOP,              # identical to imp_mCNP
         -imp_M_PB, -imp_P_PB,   # identical to imp_B_PB
         -imp_eversmok,          # identical to imp_smokstat
         -imp_previousIVF,       # identical to imp_previousIUI
         -imp_gravid)            # identical to imp_educ1

# 7. Create specific gravity-adjusted EDC variables and add study_id
# These variable names match the rename() function in 4. AFC_IF_scores_analysis.R
a_imputed <- a_imputed %>%
  mutate(
    # Add study ID as first column
    study_id = study_id,

    # Phthalate metabolites (adjusted by sgratio_pht) - 12 variables
    MBP = mBP * sgratio_pht,
    MiBP = miBP * sgratio_pht,
    MCNP = mCNP * sgratio_pht,
    MCOP = mCOP * sgratio_pht,
    MECPP = mECPP * sgratio_pht,
    MEHHP = mEHHP * sgratio_pht,
    MEHP = mEHP * sgratio_pht,
    MEOHP = mEOHP * sgratio_pht,
    MCPP = mCPP * sgratio_pht,
    MEP = mEP1 * sgratio_pht,
    MBzP = mBZP1 * sgratio_pht,
    sumDEHP = dehp * sgratio_pht,
    
    # Phenols and parabens (adjusted by sgratio_bpa) - 4 variables
    BPA = BPA * sgratio_bpa,
    BP = B_PB * sgratio_bpa,
    MP = M_PB * sgratio_bpa,
    PP = P_PB * sgratio_bpa

    # Note: BPS, BPF, BP3, TCS excluded due to >40% missing
    # Note: OPFR metabolites (BDCIPP, DPHP, ipPPP) excluded - sgratio_opfr has 74.7% missing
  ) %>%
  relocate(study_id, .before = 1)  # Ensure study_id is first column

head(a_imputed)

# 8. Save imputed dataset
imputed_data <- tibble(left_join(a_imputed, imputation_indicators, by = "study_id")) %>% 
  mutate(DOR = as.factor(as.numeric(sartnew2 == "DOR")), 
         across(c(smokstat, races, educ1), as.factor)) %>% 
  select(-B_PB, -M_PB, -P_PB, -sartnew2)

head(imputed_data)

save(imputed_data, 
     file = here("data", "imputed_EARTH.Rdata"))

# 9. Summary of imputation
cat("\n=== Imputation Summary ===\n")
cat("Original dataset dimensions:", dim(a_), "\n")
cat("Imputed dataset dimensions:", dim(a_imputed), "\n")
cat("Variables with missing data:", sum(missing_summary$n_missing > 0), "\n")
cat("Total missing values imputed:", sum(missing_summary$n_missing), "\n")
cat("NRMSE (continuous vars):", round(imputation_error$NRMSE, 4), "\n")
cat("PFC (categorical vars):", round(imputation_error$PFC, 4), "\n")