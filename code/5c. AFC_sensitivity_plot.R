pacman::p_load(
  rio,
  here,
  tidyverse
)

# ==============================================================================
# SENSITIVITY ANALYSIS PLOTTING SCRIPT
# ==============================================================================
# This program creates a test_statistic_scatter plot showing:
#   - Baseline analysis (black dots)
#   - Outliers removed sensitivity (blue dots)
#   - Extreme values imputed sensitivity (grey dots)
# Style matches the original teststat_scatter.png exactly
# ==============================================================================

cat("\n=== Creating Sensitivity Analysis Plot ===\n\n")

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

cat("Loading results...\n")

# Load baseline results from main analysis (program 4)
baseline <- readRDS(here("output", "linear_projection_output.rds")) %>%
  mutate(analysis = "Baseline")

# Load outliers removed results (program 5a)
outliers_removed <- readRDS(here("output", "sensitivity_outliers_removed.rds")) %>%
  mutate(analysis = "Outliers Removed")

# Load imputation results (program 5b)
imputation <- readRDS(here("output", "sensitivity_imputation.rds")) %>%
  mutate(analysis = "Extreme Values Imputed")

# Combine results
all_results <- bind_rows(baseline, outliers_removed, imputation)

cat("  Loaded", nrow(baseline), "baseline results\n")
cat("  Loaded", nrow(outliers_removed), "outliers removed results\n")
cat("  Loaded", nrow(imputation), "imputation results\n")

# ==============================================================================
# CREATE PLOT
# ==============================================================================

cat("\nCreating plot...\n")

# Create color mapping (baseline = black, outliers removed = blue, imputation = grey)
all_results <- all_results %>%
  mutate(
    color_group = factor(analysis, levels = c("Baseline", "Outliers Removed", "Extreme Values Imputed"))
  )

# Create the plot matching teststat_scatter.png style exactly
p <- ggplot(all_results, aes(x = statistic, y = reorder(term, statistic))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(aes(shape = Type, color = color_group), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = c("Baseline" = "black",
                                 "Outliers Removed" = "#2E86AB",
                                 "Extreme Values Imputed" = "#808080")) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0)) +
  labs(color = NULL, shape = NULL)

# Display the plot
print(p)

cat("\nPlot generated successfully.\n")
cat("\nTo save the plot, use:\n")
cat("  ggsave(filename = here('figures', 'teststat_scatter_sensitivity.png'),\n")
cat("         width = 16, height = 16, units = 'cm', dpi = 300)\n")
