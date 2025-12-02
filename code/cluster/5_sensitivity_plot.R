userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)
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

cat("\n=== Creating Sensitivity Analysis Plot ===\n\n")

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
# CREATE CONSISTENT Y-AXIS ORDERING
# ==============================================================================

# Get chemical ordering from baseline unconditional test statistics
# This will be used for all plots to maintain consistency
baseline_order <- baseline %>%
  filter(Type == "unconditional") %>%
  arrange(statistic) %>%
  pull(term)

cat("Chemical ordering based on baseline unconditional test statistics:\n")
cat(paste(baseline_order, collapse = ", "), "\n\n")

# ==============================================================================
# COMBINED SENSITIVITY PLOT - ALL SCENARIOS (UNCONDITIONAL)
# ==============================================================================

cat("Creating combined sensitivity plot for all scenarios (unconditional)...\n")

# Combine all scenarios - baseline plus all sensitivity scenarios
# Filter to unconditional models only
all_scenarios <- bind_rows(
  baseline %>% filter(Type == "unconditional"),
  sensitivity %>% filter(Type == "unconditional")
) %>%
  # Create descriptive scenario labels
  mutate(
    scenario_label = case_when(
      scenario_id == 0 ~ "Baseline (Age≥35)",
      TRUE ~ paste0(
        "Age≥", age_threshold,
        ifelse(afc_impute, ", AFC>30 Imp", ""),
        ifelse(edc_impute, ", EDC>97.5th Imp", "")
      )
    )
  ) %>%
  # Order scenarios: baseline first, then by age threshold, then imputation
  mutate(
    scenario_order = case_when(
      scenario_id == 0 ~ 0,
      TRUE ~ age_threshold * 100 + afc_impute * 10 + edc_impute * 1
    )
  ) %>%
  arrange(scenario_order) %>%
  # Apply consistent chemical ordering
  mutate(term = factor(term, levels = baseline_order))

# Get total number of scenarios
n_scenarios <- n_distinct(all_scenarios$scenario_label)
cat("  Total scenarios to plot:", n_scenarios, "\n")

# Create ordered factor for legend
scenario_levels <- all_scenarios %>%
  distinct(scenario_label, scenario_order) %>%
  arrange(scenario_order) %>%
  pull(scenario_label)

all_scenarios <- all_scenarios %>%
  mutate(scenario_label = factor(scenario_label, levels = scenario_levels))

# Create plot with all scenarios
# Use term (now a factor with baseline_order levels) directly without reorder()
p_combined <- ggplot(all_scenarios, aes(x = statistic, y = term,
                                         color = scenario_label)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_viridis_d(option = "turbo") +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

ggsave(filename = here("figures", "sensitivity_combined_unconditional.png"),
       plot = p_combined, width = 16, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# TORNADO DIAGRAM - QUANTIFY SENSITIVITY BY DIMENSION
# ==============================================================================

cat("\nCreating tornado diagram to quantify sensitivity by dimension...\n")

# Calculate changes from baseline for each dimension
# Get baseline statistics for comparison
baseline_stats <- baseline %>%
  filter(Type == "unconditional") %>%
  select(term, baseline_stat = statistic)

# AFC imputation effect (age=35, varying AFC imputation, EDC at baseline)
afc_effects <- sensitivity %>%
  filter(Type == "unconditional", age_threshold == 35, afc_impute == TRUE, edc_impute == FALSE) %>%
  select(term, statistic) %>%
  left_join(baseline_stats, by = "term") %>%
  mutate(
    change = statistic - baseline_stat,
    dimension = "AFC>30 Imputed",
    dimension_type = "AFC Imputation",
    dimension_order = 1
  )

# EDC imputation effect (age=35, AFC at baseline, varying EDC imputation)
edc_effects <- sensitivity %>%
  filter(Type == "unconditional", age_threshold == 35, afc_impute == FALSE, edc_impute == TRUE) %>%
  select(term, statistic) %>%
  left_join(baseline_stats, by = "term") %>%
  mutate(
    change = statistic - baseline_stat,
    dimension = "EDC>97.5th Imputed",
    dimension_type = "EDC Imputation",
    dimension_order = 2
  )

# Combine all effects (removed age_effects)
tornado_data <- bind_rows(afc_effects, edc_effects) %>%
  mutate(
    dimension_type = factor(dimension_type,
                           levels = c("AFC Imputation", "EDC Imputation"))
  ) %>%
  # Apply consistent chemical ordering (same as other plots)
  mutate(term = factor(term, levels = baseline_order))

# Create tornado plot (black and white)
p_tornado <- ggplot(tornado_data, aes(x = change, y = term, fill = dimension_type)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("AFC Imputation" = "gray50",
                               "EDC Imputation" = "gray80")) +
  xlab("Change in Test Statistic from Baseline") +
  ylab("Environmental Chemical") +
  labs(fill = "Sensitivity\nDimension") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

ggsave(filename = here("figures", "sensitivity_tornado.png"),
       plot = p_tornado, width = 20, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# COMBINED PLOT: BASELINE + SENSITIVITY RANGES
# ==============================================================================

cat("\nCreating combined baseline + sensitivity ranges plot...\n")

# Calculate ranges for each dimension
# AFC imputation range (comparing with/without AFC imputation at age=35, no EDC)
afc_range <- bind_rows(
  baseline %>% filter(Type == "unconditional") %>% select(term, statistic),
  sensitivity %>% filter(Type == "unconditional", age_threshold == 35, afc_impute == TRUE, edc_impute == FALSE) %>%
    select(term, statistic)
) %>%
  group_by(term) %>%
  summarize(
    min_stat = min(statistic),
    max_stat = max(statistic),
    .groups = "drop"
  ) %>%
  left_join(baseline_stats, by = "term") %>%
  mutate(
    dimension = "AFC Imputation",
    lower = min_stat,
    upper = max_stat,
    point = baseline_stat
  )

# EDC imputation range (comparing with/without EDC imputation at age=35, no AFC)
edc_range <- bind_rows(
  baseline %>% filter(Type == "unconditional") %>% select(term, statistic),
  sensitivity %>% filter(Type == "unconditional", age_threshold == 35, afc_impute == FALSE, edc_impute == TRUE) %>%
    select(term, statistic)
) %>%
  group_by(term) %>%
  summarize(
    min_stat = min(statistic),
    max_stat = max(statistic),
    .groups = "drop"
  ) %>%
  left_join(baseline_stats, by = "term") %>%
  mutate(
    dimension = "EDC Imputation",
    lower = min_stat,
    upper = max_stat,
    point = baseline_stat
  )

# Combine all ranges (removed age_range)
range_data <- bind_rows(afc_range, edc_range) %>%
  mutate(
    dimension = factor(dimension, levels = c("AFC Imputation", "EDC Imputation")),
    term = factor(term, levels = baseline_order)
  )

# ==============================================================================
# CREATE SHARED SPACING FOR BOTH PLOTS (UNCONDITIONAL AND CONDITIONAL)
# ==============================================================================
# This ensures both plots have identical y-axis spacing
spacing_factor <- 2.5
# Calculate y-axis breaks (one per chemical, in baseline order)
# This is used for BOTH unconditional and conditional plots
y_breaks <- (seq_along(baseline_order) - 1) * spacing_factor + 1

# Create plot with baseline points and sensitivity ranges
# Use baseline_order directly in the aesthetic mapping to ensure consistent spacing
p_combined_ranges <- ggplot(range_data, aes(x = point, group = dimension)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_errorbar(aes(xmin = lower, xmax = upper,
                    y = ((as.numeric(term) - 1) * spacing_factor + 1) + (as.numeric(dimension) - 1.5) * 0.65),
                width = 0, linewidth = 0.5, color = "black", orientation = "y") +
  geom_point(aes(x = point,
                 y = ((as.numeric(term) - 1) * spacing_factor + 1) + (as.numeric(dimension) - 1.5) * 0.65,
                 shape = dimension),
             size = 2.5, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c("AFC Imputation" = 22,      # square
                                "EDC Imputation" = 24)) +   # triangle
  scale_y_continuous(breaks = y_breaks, labels = baseline_order) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(shape = "Sensitivity\nDimension") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

ggsave(filename = here("figures", "sensitivity_baseline_ranges_unconditional.png"),
       plot = p_combined_ranges, width = 20, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# CONDITIONAL VERSION: BASELINE + SENSITIVITY RANGES (CONDITIONAL)
# ==============================================================================

cat("\nCreating conditional baseline + sensitivity ranges plot...\n")

# Get baseline conditional statistics for comparison
baseline_stats_cond <- baseline %>%
  filter(Type == "conditional") %>%
  select(term, baseline_stat = statistic)

# Calculate ranges for each dimension - CONDITIONAL
# AFC imputation range (comparing with/without AFC imputation at age=35, no EDC)
afc_range_cond <- bind_rows(
  baseline %>% filter(Type == "conditional") %>% select(term, statistic),
  sensitivity %>% filter(Type == "conditional", age_threshold == 35, afc_impute == TRUE, edc_impute == FALSE) %>%
    select(term, statistic)
) %>%
  group_by(term) %>%
  summarize(
    min_stat = min(statistic),
    max_stat = max(statistic),
    .groups = "drop"
  ) %>%
  left_join(baseline_stats_cond, by = "term") %>%
  mutate(
    dimension = "AFC Imputation",
    lower = min_stat,
    upper = max_stat,
    point = baseline_stat
  )

# EDC imputation range (comparing with/without EDC imputation at age=35, no AFC)
edc_range_cond <- bind_rows(
  baseline %>% filter(Type == "conditional") %>% select(term, statistic),
  sensitivity %>% filter(Type == "conditional", age_threshold == 35, afc_impute == FALSE, edc_impute == TRUE) %>%
    select(term, statistic)
) %>%
  group_by(term) %>%
  summarize(
    min_stat = min(statistic),
    max_stat = max(statistic),
    .groups = "drop"
  ) %>%
  left_join(baseline_stats_cond, by = "term") %>%
  mutate(
    dimension = "EDC Imputation",
    lower = min_stat,
    upper = max_stat,
    point = baseline_stat
  )

# Combine all ranges (removed age_range_cond)
range_data_cond <- bind_rows(afc_range_cond, edc_range_cond) %>%
  mutate(
    dimension = factor(dimension, levels = c("AFC Imputation", "EDC Imputation")),
    term = factor(term, levels = baseline_order)
  )

# Create conditional plot (identical format to unconditional)
# IMPORTANT: Uses same spacing_factor and y_breaks defined earlier
# Calculate y-position inline using the same formula as unconditional plot
p_combined_ranges_cond <- ggplot(range_data_cond, aes(x = point, group = dimension)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_errorbar(aes(xmin = lower, xmax = upper,
                    y = ((as.numeric(term) - 1) * spacing_factor + 1) + (as.numeric(dimension) - 1.5) * 0.65),
                width = 0, linewidth = 0.5, color = "black", orientation = "y") +
  geom_point(aes(x = point,
                 y = ((as.numeric(term) - 1) * spacing_factor + 1) + (as.numeric(dimension) - 1.5) * 0.65,
                 shape = dimension),
             size = 2.5, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c("AFC Imputation" = 22,      # square
                                "EDC Imputation" = 24)) +   # triangle
  scale_y_continuous(breaks = y_breaks, labels = baseline_order) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  labs(shape = "Sensitivity\nDimension") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

ggsave(filename = here("figures", "sensitivity_baseline_ranges_conditional.png"),
       plot = p_combined_ranges_cond, width = 20, height = 16, units = "cm", dpi = 300)

# ==============================================================================
# COMBINED TWO-PANEL FIGURE: UNCONDITIONAL (A) AND CONDITIONAL (B)
# ==============================================================================

cat("\nCreating combined two-panel figure (unconditional + conditional)...\n")

# Modify plots to remove legends (will create shared legend)
p_combined_ranges_A <- p_combined_ranges +
  labs(title = "A. Unconditional Model") +
  theme(legend.position = "none")

p_combined_ranges_B <- p_combined_ranges_cond +
  ylab(NULL) +  # Remove y-axis label from panel B
  labs(title = "B. Conditional Model") +
  theme(legend.position = "none")

# Extract legend from one of the plots (they're identical)
library(gridExtra)
library(grid)

# Get legend from the original plot
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

shared_legend <- get_legend(p_combined_ranges)

# Combine into two-panel figure with shared legend
p_two_panel <- grid.arrange(
  arrangeGrob(p_combined_ranges_A, p_combined_ranges_B, ncol = 2),
  shared_legend,
  ncol = 2,
  widths = c(10, 1)  # Main plots take 10 units, legend takes 1 unit
)

ggsave(filename = here("figures", "sensitivity_baseline_ranges_combined.png"),
       plot = p_two_panel, width = 42, height = 16, units = "cm", dpi = 300)

cat("  Two-panel figure saved: figures/sensitivity_baseline_ranges_combined.png\n")

# Calculate and save summary statistics
sensitivity_summary <- tornado_data %>%
  group_by(dimension_type) %>%
  summarize(
    mean_abs_change = mean(abs(change)),
    median_abs_change = median(abs(change)),
    max_abs_change = max(abs(change)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_change))

saveRDS(sensitivity_summary, file = here("output", "sensitivity_summary_by_dimension.rds"))

cat("\n=== Sensitivity plots created successfully ===\n")
cat("Generated figures:\n")
cat("  1. figures/sensitivity_combined_unconditional.png - All scenarios (no legend)\n")
cat("  2. figures/sensitivity_tornado.png - Changes by dimension (B&W)\n")
cat("  3. figures/sensitivity_baseline_ranges_unconditional.png - Unconditional + ranges\n")
cat("  4. figures/sensitivity_baseline_ranges_conditional.png - Conditional + ranges\n")
cat("  5. figures/sensitivity_baseline_ranges_combined.png - Two-panel (A & B)\n")
cat("\nSummary statistics by dimension:\n")
print(sensitivity_summary)
cat("\n")
