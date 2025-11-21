pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  reshape2,
  dbscan,      # For LOF (Local Outlier Factor)
  factoextra,  # For PCA visualization
  MVN          # For multivariate normality and outlier tests
)

## goal of this program is to explore distribution of environmental exposures 
## in EARTH data
## 
## Compare both trimmed and untrimmed versions

load(here("data", "afc_clean_trunc.Rdata"))
load(here("data", "afc_clean_notrunc.Rdata"))

## some basic outlier / leverage analyses
## do the people with extreme chem values also have large outcome values?
##
# 17 EDC variables with <40% missing (16 SG-adjusted + 1 Hg)
# (defined in file 1. AFC_data_man.R)
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "sumDEHP",
              "BPA", "BP", "MP", "PP",
              "Hg")

plot_distributions <- function(data, variables, 
                               ncol = 4, title = NULL,
                               linewidth = .5) {
  # Create long format data
  plot_data_long <- data %>%
    select(all_of(variables)) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    drop_na()
  
  # Create faceted density + histogram plot
  p <- ggplot(plot_data_long) +
    geom_histogram(aes(x = value,
                       y = after_stat(density)),
                   bins = 30,
                   fill = "white",
                   color = "black",
                   linewidth = linewidth) +
    geom_density(aes(x = value),
                 kernel = "epanechnikov",
                 bw = "SJ",
                 adjust = 1,
                 linewidth = linewidth,
                 color = "#2297E6") +
    facet_wrap(~variable, scales = "free", ncol = ncol) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("Density") +
    xlab("Value") +
    theme_classic() +
    theme(strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold"))
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

plot_vars <- c("AFCt", "age", env_vars)

p_notrunc <- plot_distributions(
  data = afc_clean_notrunc,
  variables = plot_vars,
  ncol = 4,
  title = "Distribution of AFC, age, and Environmental Chemicals (Non-truncated)"
)
print(p_notrunc)

ggsave(here("figures", "env_vars_notrunc.png"), units = "cm", width = 16, height = 16, dpi = 300)

p_trunc <- plot_distributions(
  data = afc_clean_trunc,
  variables = plot_vars,
  ncol = 4,
  title = "Distribution of AFC, age, and Environmental Chemicals (Truncated)"
)
print(p_trunc)

ggsave(here("figures", "env_vars_trunc.png"), units = "cm", width = 16, height = 16, dpi = 300)

## multivariate outlier detection
# Prepare multivariate data (outcome + all environmental variables)
multivar_data <- afc_clean_notrunc %>%
  select(AFCt, age, all_of(env_vars)) 

cat("Sample size for outlier detection:", nrow(multivar_data), "\n")
cat("Variables included:", ncol(multivar_data), "\n\n")

# Store row indices to match back to original data
multivar_data$row_id <- as.numeric(rownames(multivar_data))

# ============================================================================
# METHOD 1: Mahalanobis Distance
# ============================================================================
# Measures how far each observation is from the center of the multivariate
# distribution, accounting for correlations between variables
# $$D_M(\mathbf{x}) = \sqrt{(\mathbf{x} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\mu})}$$


cat("--- Method 1: Mahalanobis Distance ---\n")

# Calculate Mahalanobis distance
mahala_dist <- mahalanobis(
  x = multivar_data[, c("AFCt", "age", env_vars)],
  center = colMeans(multivar_data[, c("AFCt", "age", env_vars)]),
  cov = cov(multivar_data[, c("AFCt", "age", env_vars)])
)

# Chi-square critical value (p < 0.001 for conservative threshold)
# df = number of variables
crit_value <- qchisq(0.999, df = length(c("AFCt", "age", env_vars)))

multivar_data$mahala_dist <- mahala_dist
multivar_data$mahala_outlier <- mahala_dist > crit_value

cat("Mahalanobis outliers (p < 0.001):", sum(multivar_data$mahala_outlier), "\n")
cat("Critical value:", round(crit_value, 2), "\n\n")

plot_data <- tibble(mahala_dist)

# Plot 1: Mahalanobis distance
p1 <- ggplot(multivar_data, aes(x = 1:nrow(multivar_data), y = mahala_dist)) +
  geom_point(aes(color = mahala_outlier), alpha = 0.6) +
  geom_hline(yintercept = crit_value, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray50", "red")) +
  labs(title = "Mahalanobis Distance",
       x = "Observation",
       y = "Mahalanobis Distance",
       color = "Outlier") +
  theme_classic()
print(p1)

# ============================================================================
# METHOD 2: Principal Component Analysis (PCA) Based Outliers
# ============================================================================
# Identifies outliers in the reduced dimensional space

cat("--- Method 2: PCA-Based Outliers ---\n")

# Perform PCA (scale variables for fair comparison)
pca_result <- prcomp(multivar_data[, c("AFCt", "age", env_vars)],
                     scale. = TRUE,
                     center = TRUE)

cat("Variance explained by first 5 PCs:",
    round(sum(summary(pca_result)$importance[2, 1:5]) * 100, 1), "%\n")


# Extract variance explained
pca_var <- summary(pca_result)$importance
pca_df <- data.frame(
  PC = paste0("PC", 1:ncol(pca_var)),
  PC_num = 1:ncol(pca_var),
  variance = pca_var[2, ] * 100,  # Proportion of variance as percentage
  cumulative = pca_var[3, ] * 100  # Cumulative proportion as percentage
)

# Scree plot with cumulative variance
ggplot(pca_df) +
  geom_line(aes(x = PC_num, y = variance, color = "Individual"),
            linewidth = 1) +
  geom_point(aes(x = PC_num, y = variance, color = "Individual"),
             size = 3) +
  geom_line(aes(x = PC_num, y = cumulative, color = "Cumulative"),
            linewidth = 1) +
  geom_point(aes(x = PC_num, y = cumulative, color = "Cumulative"),
             size = 3) +
  geom_hline(yintercept = 80,
             linetype = "dashed",
             color = "gray50",
             linewidth = 0.5) +
  scale_color_manual(values = c("Individual" = "#2297E6",
                                "Cumulative" = "#DF536B")) +
  scale_x_continuous(breaks = 1:min(10, ncol(pca_var))) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 105)) +
  labs(title = "PCA Variance Explained",
       x = "Principal Component",
       y = "Variance Explained (%)",
       color = "Type") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# Get PC scores
pc_scores <- as.data.frame(pca_result$x[, 1:5])

# Calculate distance in PC space (using first 5 PCs)
pc_dist <- sqrt(rowSums(pc_scores^2))
pc_threshold <- quantile(pc_dist, 0.99)  # Top 1% as outliers

multivar_data$pc_dist <- pc_dist
multivar_data$pc_outlier <- pc_dist > pc_threshold

cat("PCA outliers (top 1%):", sum(multivar_data$pc_outlier), "\n")
cat("PC distance threshold:", round(pc_threshold, 2), "\n\n")

# Plot 2: PCA biplot with outliers highlighted
p2 <- ggplot(data.frame(PC1 = pca_result$x[,1],
                        PC2 = pca_result$x[,2],
                        outlier = multivar_data$pc_outlier)) +
  geom_point(aes(x = PC1, y = PC2, color = outlier), alpha = 0.6) +
  scale_color_manual(values = c("gray50", "red")) +
  labs(title = "PCA: First Two Components",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
       color = "Outlier") +
  theme_classic()
plot(p2)

# ============================================================================
# METHOD 3: Local Outlier Factor (LOF)
# ============================================================================
# Density-based method that identifies observations in low-density regions

cat("--- Method 3: Local Outlier Factor (LOF) ---\n")

# Calculate LOF scores (k = 20 neighbors)
lof_scores <- lof(multivar_data[, c("AFCt", env_vars)], k = 20)

# LOF > 1.5 suggests outlier (observations in sparser regions)
lof_threshold <- 1.5

multivar_data$lof_score <- lof_scores
multivar_data$lof_outlier <- lof_scores > lof_threshold

cat("LOF outliers (LOF > 1.5):", sum(multivar_data$lof_outlier), "\n")
cat("Mean LOF score:", round(mean(lof_scores), 3), "\n")
cat("Max LOF score:", round(max(lof_scores), 3), "\n\n")

# Plot 3: LOF scores
p3 <- ggplot(multivar_data, aes(x = 1:nrow(multivar_data), y = lof_score)) +
  geom_point(aes(color = lof_outlier), alpha = 0.6) +
  geom_hline(yintercept = lof_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray50", "red")) +
  scale_y_log10() +
  labs(title = "Local Outlier Factor (LOF)",
       x = "Observation",
       y = "LOF Score (log scale)",
       color = "Outlier") +
  theme_classic()

# ============================================================================
# METHOD 4: Regression Diagnostics 
# ============================================================================
# Cook's distance and leverage from regression model

cat("--- Method 4: Regression Diagnostics ---\n")

# Create formula: AFC_t ~ all chemicals
formula_vars <- paste(env_vars, collapse = " + ")
model_formula <- as.formula(paste0("AFCt ~ ", formula_vars))

# Fit linear model
mod <- lm(model_formula, data = afc_clean_notrunc)

# Get diagnostic data with broom
diagnostic_data <- broom::augment(mod)

# Merge regression diagnostics
reg_diag <- diagnostic_data %>%
  mutate(row_id = as.numeric(rownames(diagnostic_data))) %>%
  select(row_id, .cooksd, .hat, .std.resid)

multivar_data <- multivar_data %>%
  left_join(reg_diag, by = "row_id")

# Flag high Cook's D (> 4/n) and high leverage (> 2*p/n)
n <- nrow(multivar_data)
p <- length(env_vars) + 1
cooksd_threshold <- 4 / n
leverage_threshold <- 2 * p / n

multivar_data$high_cooksd <- multivar_data$.cooksd > cooksd_threshold
multivar_data$high_leverage <- multivar_data$.hat > leverage_threshold
multivar_data$high_stdresid <- abs(multivar_data$.std.resid) > 3

cat("High Cook's D (> 4/n):", sum(multivar_data$high_cooksd, na.rm = TRUE), "\n")
cat("High leverage (> 2p/n):", sum(multivar_data$high_leverage, na.rm = TRUE), "\n")
cat("High std residuals (|z| > 3):", sum(multivar_data$high_stdresid, na.rm = TRUE), "\n\n")

p4 <- ggplot(multivar_data, aes(x = 1:nrow(multivar_data), y = .cooksd)) +
  geom_point(aes(color = high_cooksd), alpha = 0.6) +
  geom_hline(yintercept = cooksd_threshold,
             linetype = "dashed",
             color = "red") +
  scale_color_manual(values = c("gray50", "red"),
                     labels = c("Normal", "High Influence")) +
  labs(title = "Cook's Distance",
       x = "Observation",
       y = "Cook's Distance",
       color = "Influence") +
  theme_classic()
print(p4)

# ============================================================================
# SUMMARY: Flagged by Multiple Methods
# ============================================================================

cat("=== OUTLIER SUMMARY ===\n")

# Count how many methods flagged each observation
multivar_data$n_flags <- rowSums(multivar_data[, c("mahala_outlier",
                                                     "pc_outlier",
                                                     "lof_outlier",
                                                     "high_cooksd")])

# Identify observations flagged by 2+ methods
consensus_outliers <- multivar_data %>%
  filter(n_flags >= 2) %>%
  arrange(desc(n_flags))

cat("Observations flagged by 2+ methods:", nrow(consensus_outliers), "\n")
cat("Observations flagged by 3+ methods:", sum(multivar_data$n_flags >= 3), "\n")
cat("Observations flagged by all 4 methods:", sum(multivar_data$n_flags == 4), "\n\n")

# ============================================================================
# VISUALIZATIONS
# ============================================================================

# Combine plots
outlider_grid <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(plot = outlider_grid, 
       here("figures", "outlier_grid.png"), units = "cm", width = 16, height = 16, dpi = 300)

# Plot 5: Heatmap of flagged observations
if (nrow(consensus_outliers) > 0) {
  consensus_matrix <- consensus_outliers %>%
    select(row_id, mahala_outlier, pc_outlier, lof_outlier, high_cooksd) %>%
    pivot_longer(cols = -row_id, names_to = "method", values_to = "flagged")

  p5 <- ggplot(consensus_matrix, aes(x = method, y = factor(row_id), fill = flagged)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("white", "red")) +
    labs(title = "Consensus Outliers (Flagged by 2+ Methods)",
         x = "Detection Method",
         y = "Observation ID",
         fill = "Flagged") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p5)
}

ggsave(here("figures", "outlier_heatmap.png"), units = "cm", width = 16, height = 16, dpi = 300)

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save outlier analysis results
outlier_results <- multivar_data %>%
  select(row_id, AFCt, age, mahala_dist, mahala_outlier,
         pc_dist, pc_outlier, lof_score, lof_outlier,
         .cooksd, high_cooksd, n_flags) %>%
  arrange(desc(n_flags), desc(mahala_dist))

# Save to misc directory for further investigation
write_csv(outlier_results, here("misc", "outlier_analysis_results.csv"))

cat("\n=== Results saved to misc/outlier_analysis_results.csv ===\n")

# Display top consensus outliers
cat("\nTop 10 observations flagged by multiple methods:\n")
outliers_id <- outlier_results %>%
        filter(n_flags >= 2) %>%
        head(10) %>%
        select(row_id, AFCt, age, n_flags, mahala_outlier, pc_outlier, lof_outlier, high_cooksd)


saveRDS(outliers_id, here("output", "outliers_id.rds"))

