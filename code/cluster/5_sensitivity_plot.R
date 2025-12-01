pacman::p_load(
  rio,
  here,
  tidyverse
)

# ==============================================================================
# SENSITIVITY ANALYSIS PLOTTING SCRIPT
# ==============================================================================
# This program creates visualization of sensitivity analysis results across:
#   - Age thresholds (35-40 years)
#   - AFC>30 imputation (no vs yes)
#   - EDC>97.5th percentile imputation (no vs yes)
# ==============================================================================

cat("\n=== Creating Sensitivity Analysis Plots ===\n\n")

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

cat("Loading results...\n")

# Load baseline results from main analysis (program 4)
baseline <- readRDS(here("output", "linear_projection_output.rds")) %>%
  mutate(
    age_threshold = 35,
    afc_impute = FALSE,
    edc_impute = FALSE,
    scenario_id = 0
  )

# Load sensitivity analysis results (program 5)
sensitivity <- readRDS(here("output", "sensitivity_analysis_results.rds"))

cat("  Loaded", nrow(baseline), "baseline results\n")
cat("  Loaded", nrow(sensitivity), "sensitivity results\n")
cat("  Unique scenarios:", n_distinct(sensitivity$scenario_id), "\n\n")

# ==============================================================================
# PLOT 1: AGE THRESHOLD SENSITIVITY (UNCONDITIONAL)
# ==============================================================================

cat("Creating age threshold sensitivity plot (unconditional)...\n")

# Filter to baseline imputation settings (no AFC/EDC imputation) and unconditional
age_sensitivity <- sensitivity %>%
  filter(
    afc_impute == FALSE,
    edc_impute == FALSE,
    Type == "unconditional"
  ) %>%
  bind_rows(baseline %>% filter(Type == "unconditional"))

p1 <- ggplot(age_sensitivity, aes(x = statistic, y = reorder(term, statistic), color = factor(age_threshold))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(color = "Age Threshold") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

ggsave(filename = here("figures", "sensitivity_age_threshold.png"),
       plot = p1, width = 20, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# PLOT 2: AFC IMPUTATION SENSITIVITY (UNCONDITIONAL)
# ==============================================================================

cat("Creating AFC imputation sensitivity plot (unconditional)...\n")

# Filter to age=35, no EDC imputation, unconditional
afc_sensitivity <- sensitivity %>%
  filter(
    age_threshold == 35,
    edc_impute == FALSE,
    Type == "unconditional"
  ) %>%
  bind_rows(baseline %>% filter(Type == "unconditional")) %>%
  mutate(afc_impute_label = ifelse(afc_impute, "AFC>30 Imputed", "Baseline"))

p2 <- ggplot(afc_sensitivity, aes(x = statistic, y = reorder(term, statistic), color = afc_impute_label)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(size = 3) +
  scale_color_manual(values = c("Baseline" = "black", "AFC>30 Imputed" = "#E74C3C")) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(color = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0))

ggsave(filename = here("figures", "sensitivity_afc_imputation.png"),
       plot = p2, width = 16, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# PLOT 3: EDC IMPUTATION SENSITIVITY (UNCONDITIONAL)
# ==============================================================================

cat("Creating EDC imputation sensitivity plot (unconditional)...\n")

# Filter to age=35, no AFC imputation, unconditional
edc_sensitivity <- sensitivity %>%
  filter(
    age_threshold == 35,
    afc_impute == FALSE,
    Type == "unconditional"
  ) %>%
  bind_rows(baseline %>% filter(Type == "unconditional")) %>%
  mutate(edc_impute_label = ifelse(edc_impute, "EDC>97.5th Imputed", "Baseline"))

p3 <- ggplot(edc_sensitivity, aes(x = statistic, y = reorder(term, statistic), color = edc_impute_label)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(size = 3) +
  scale_color_manual(values = c("Baseline" = "black", "EDC>97.5th Imputed" = "#3498DB")) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(color = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0))

ggsave(filename = here("figures", "sensitivity_edc_imputation.png"),
       plot = p3, width = 16, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# PLOT 4: COMBINED SENSITIVITY (UNCONDITIONAL, SELECTED SCENARIOS)
# ==============================================================================

cat("Creating combined sensitivity plot (unconditional)...\n")

# Select key scenarios to compare
combined_sensitivity <- bind_rows(
  baseline %>% filter(Type == "unconditional") %>% mutate(scenario_label = "Baseline (Age≥35)"),
  sensitivity %>% filter(age_threshold == 35, afc_impute == FALSE, edc_impute == TRUE, Type == "unconditional") %>%
    mutate(scenario_label = "EDC>97.5th Imputed"),
  sensitivity %>% filter(age_threshold == 35, afc_impute == TRUE, edc_impute == FALSE, Type == "unconditional") %>%
    mutate(scenario_label = "AFC>30 Imputed"),
  sensitivity %>% filter(age_threshold == 40, afc_impute == FALSE, edc_impute == FALSE, Type == "unconditional") %>%
    mutate(scenario_label = "Age≥40"),
  sensitivity %>% filter(age_threshold == 35, afc_impute == TRUE, edc_impute == TRUE, Type == "unconditional") %>%
    mutate(scenario_label = "Both Imputed")
)

p4 <- ggplot(combined_sensitivity, aes(x = statistic, y = reorder(term, statistic),
                                        color = factor(scenario_label, levels = c("Baseline (Age≥35)", "Age≥40", "AFC>30 Imputed", "EDC>97.5th Imputed", "Both Imputed")))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Baseline (Age≥35)" = "black",
                                "Age≥40" = "#9B59B6",
                                "AFC>30 Imputed" = "#E74C3C",
                                "EDC>97.5th Imputed" = "#3498DB",
                                "Both Imputed" = "#2ECC71")) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(color = "Scenario") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

ggsave(filename = here("figures", "sensitivity_combined.png"),
       plot = p4, width = 24, height = 16, units = "cm", dpi = 300)

cat("\n=== All sensitivity plots created successfully ===\n")
cat("Generated 4 plots:\n")
cat("  1. sensitivity_age_threshold.png\n")
cat("  2. sensitivity_afc_imputation.png\n")
cat("  3. sensitivity_edc_imputation.png\n")
cat("  4. sensitivity_combined.png\n")
