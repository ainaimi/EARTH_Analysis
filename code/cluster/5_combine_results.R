#!/usr/bin/env Rscript
userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)
# ==============================================================================
# COMBINE SENSITIVITY ANALYSIS RESULTS FROM CLUSTER
# ==============================================================================
# This script combines the output from all cluster jobs into a single file.
# Run this AFTER all cluster jobs have completed successfully.
#
# Usage:
#   Rscript code/5_combine_results.R
# ==============================================================================

pacman::p_load(rio, here, tidyverse)

cat("\n")
cat("================================================================================\n")
cat("COMBINING SENSITIVITY ANALYSIS RESULTS\n")
cat("================================================================================\n\n")

# ==============================================================================
# FIND AND LOAD ALL TASK RESULTS
# ==============================================================================

cat("Searching for task result files...\n")

# Find all task result files
task_files <- list.files(
  path = here("output"),
  pattern = "^sensitivity_task_[0-9]+\\.rds$",
  full.names = TRUE
)

n_tasks <- length(task_files)

if (n_tasks == 0) {
  stop("Error: No task result files found in output/ directory.\n",
       "Make sure all cluster jobs have completed successfully.")
}

cat("Found", n_tasks, "task result files:\n")
for (f in task_files) {
  cat("  -", basename(f), "\n")
}
cat("\n")

# ==============================================================================
# LOAD AND COMBINE RESULTS
# ==============================================================================

cat("Loading and combining results...\n")

all_results <- map_dfr(task_files, function(f) {
  cat("  Reading:", basename(f), "\n")
  readRDS(f)
})

cat("\nCombined results:\n")
cat("  Total rows:", nrow(all_results), "\n")
cat("  Unique scenarios:", n_distinct(all_results$scenario_id), "\n")
cat("  Expected scenarios: 24\n")

# Check if we have all expected scenarios
if (n_distinct(all_results$scenario_id) != 24) {
  warning("Warning: Expected 24 scenarios but found ",
          n_distinct(all_results$scenario_id),
          " scenarios.")
  cat("\nMissing scenarios:\n")
  expected <- 1:24
  found <- unique(all_results$scenario_id)
  missing <- setdiff(expected, found)
  print(missing)
}

# ==============================================================================
# SAVE COMBINED RESULTS
# ==============================================================================

cat("\nSaving combined results...\n")

output_file <- here("output", "sensitivity_analysis_results.rds")
saveRDS(all_results, file = output_file)

cat("Combined results saved to:", output_file, "\n")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("\n================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

# Summary by scenario
scenario_summary <- all_results %>%
  group_by(scenario_id, age_threshold, afc_impute, edc_impute) %>%
  summarise(
    n_tests = n(),
    n_sig_05 = sum(abs(statistic) > 1.96),
    .groups = "drop"
  ) %>%
  arrange(scenario_id)

print(scenario_summary, n = Inf)

cat("\n================================================================================\n")
cat("COMBINATION COMPLETE\n")
cat("You can now run: Rscript code/5c. AFC_sensitivity_plot.R\n")
cat("================================================================================\n\n")
