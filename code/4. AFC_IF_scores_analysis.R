pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  reshape2,
  clubSandwich,
  mgcViz,
  mgcv,
  ggh4x,
  ggrepel,
  SuperLearner,
  gridExtra
)

load(here("data", "afc_clean_notrunc_IF.Rdata"))

if_data <- afc_clean_notrunc %>%
  select(dr_scores)

if_data <- data.frame(if_data)

names(if_data) <- c("dr_scores")

# 17 EDC variables with <40% missing (16 SG-adjusted + 1 Hg)
# (defined in file 1. AFC_data_man.R)
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", #"sumDEHP",
              "BPA", "BP", "MP", "PP",
              "Hg")

# ATEs
ate_func <- function(a){

  a <- unlist(a)

  aipw_psi <- mean(a)

  aipw_se <- sd(a)/sqrt(length(a))

  aipw_ate <- c(aipw_psi, aipw_se)

  aipw_ate <- cbind(t(aipw_ate),
                    aipw_ate[1] - 1.96*aipw_ate[2],
                    aipw_ate[1] + 1.96*aipw_ate[2])

  colnames(aipw_ate) <- c("estimate", "std.err", "LCL", "UCL")

  return(aipw_ate)
}

res_ate <- ate_func(if_data[,1])

row.names(res_ate) <- c("AIPW with Stacking")

res_ate

res_dat <- data.frame(res_ate)


GGally::ggpairs(log(afc_clean_notrunc[,env_vars]))

#### DR Learner
## conditional models
lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>%
  mutate(across(all_of(env_vars), log))
model_dr_cond <- lm(dr_scores ~ ., data = lin_proj)
dr_learner_conditional <- tidy(coeftest(model_dr_cond, vcov = vcovHC(model_dr_cond, type = "HC3"))) %>%
  filter(term != "(Intercept)")
dr_learner_conditional <- dr_learner_conditional  %>%
  rename_with(~ paste0("dr_", .), .cols = -term)
# unconditional models
dr_learner_unconditional <- map_dfr(env_vars, function(var) {
  formula_str <- paste("dr_scores ~", var)
  model <- lm(as.formula(formula_str), data = lin_proj)
  tidy(coeftest(model, vcov = vcovHC(model, type = "HC3"))) %>%
    filter(term != "(Intercept)")  # Keep only the env_var coefficient
})
dr_learner_unconditional <- dr_learner_unconditional  %>%
  rename_with(~ paste0("dr_", .), .cols = -term)

table2_conditional <- dr_learner_conditional %>%
  mutate(Type = "conditional", .before=1)

table2_unconditional <- dr_learner_unconditional %>%
  mutate(Type = "unconditional", .before=1)

table2 <- rbind(table2_unconditional, table2_conditional) 

saveRDS(table2, file = here("output", "linear_projection_output.rds"))

# Create plot of DR learner test statistics

ggplot(table2, aes(x = dr_statistic, y = reorder(term, dr_statistic))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = c(-1.96, 1.96), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = c(-1.645, 1.645), color = "gray50", linetype = "dotted") +
  geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
  geom_point(aes(shape = Type), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  xlim(-3, 3) +
  xlab("Test Statistic, Linear Projection") +
  ylab("Environmental Chemical") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0))

ggsave(filename = here("figures", "teststat_scatter.png"),
       width = 16, height = 16, units = "cm", dpi = 300)


## modifier plots for all variables
var_names <- unique(table2$term)

# Define x-axis labels for 17 environmental variables
x_labels <- c(
  "MBP" = "log(MBP), μg/L",
  "MiBP" = "log(MiBP), μg/L",
  "MCNP" = "log(MCNP), μg/L",
  "MCOP" = "log(MCOP), μg/L",
  "MECPP" = "log(MECPP), μg/L",
  "MEHHP" = "log(MEHHP), μg/L",
  "MEHP" = "log(MEHP), μg/L",
  "MEOHP" = "log(MEOHP), μg/L",
  "MCPP" = "log(MCPP), μg/L",
  "MEP" = "log(MEP), μg/L",
  "MBzP" = "log(MBzP), μg/L",
  #"sumDEHP" = "log(∑DEHP), μg/L",
  "BPA" = "log(BPA), μg/L",
  "BP" = "log(BP), μg/L",
  "MP" = "log(MP), μg/L",
  "PP" = "log(PP), μg/L",
  "Hg" = "log(Hg), ppm"
)

## linear
# Function with linear models
create_cate_plot <- function(var_name, data, x_label, ate_estimates, show_ylab = TRUE) {
  # Select and clean data for this variable
  plot_data <- data %>%
    select(dr_scores, all_of(var_name)) %>%
    mutate(across(all_of(var_name), ~log(. + 0.01)))

  # Build formula dynamically for linear models
  formula_dr_str <- paste0("dr_scores ~ ", var_name)

  # Fit linear model
  model_dr <- lm(as.formula(formula_dr_str), data = plot_data)

  # Add predictions to data
  plot_data$pred_dr <- predict(model_dr)
  plot_data$log_var <- plot_data[[var_name]]

  # Create plot
  p <- ggplot(plot_data, aes(x = log_var)) +
    geom_line(aes(y = pred_dr), color = "black") +
    geom_rug(sides = "b", length = unit(2, "mm")) +
    geom_hline(yintercept = ate_estimates, color = "gray", linetype = "dashed") +
    xlab(x_label) + ylim(-8, -2) +
    theme_classic(base_size = 10)

  # Conditionally add y-axis label
  if (show_ylab) {
    p <- p + ylab("Age-AFC Association")
  } else {
    p <- p + ylab(NULL)
  }

  return(p)
}

# Create all plots for all EDCs
var_names <- unique(table2$term)
ncol_grid <- 4
plot_list <- lapply(seq_along(var_names), function(i) {
  var <- var_names[i]
  show_ylab <- (i - 1) %% ncol_grid == 0
  create_cate_plot(var, afc_clean_notrunc, x_labels[var], res_dat$estimate,
                   show_ylab = show_ylab)
})

# Combine plots
plot_grid <- grid.arrange(grobs = plot_list, ncol = ncol_grid)

ggsave(plot = plot_grid, filename = here("figures", "cate_functions_full_linear.png"),
       width = 20, height = 20, units = "cm", dpi = 300)

## SuperLearner-based CATE function estimation
# Create custom learner wrapper for EARTH
loess_learner <- create.Learner("SL.loess", tune = list(span = seq(1,2,by=.05)))
gam_learner <- create.Learner("SL.gam", tune = list(deg.gam = seq(1,3)))
set.seed(123)
# Function with SuperLearner that uses different learners based on chemical
create_cate_plot_sl <- function(var_name, data, x_label, ate_estimates, show_ylab = TRUE) {
  # Select and clean data for this variable
  plot_data <- data %>%
    select(dr_scores, all_of(var_name)) %>%
    mutate(across(all_of(var_name), ~log(. + 0.01))) %>%
    drop_na()  # Remove any NA values

  # Extract predictor variable and sort by it
  x_vals <- plot_data[[var_name]]
  sort_order <- order(x_vals)
  x_sorted <- x_vals[sort_order]
  y_dr_sorted <- plot_data$dr_scores[sort_order]

  # Prepare data for SuperLearner
  X_train <- data.frame(x = x_sorted)

  # Select SuperLearner library based on chemical
  if (var_name %in% c("MCOP", "MEHP")) {
    SL.library <- loess_learner$names
    cat("\nUsing LOESS learner for", var_name, "\n")
  } else {
    SL.library <- gam_learner$names
    cat("\nUsing GAM learner for", var_name, "\n")
  } 

  num.folds <- 10
  folds <- sort(seq(nrow(plot_data)) %% num.folds) + 1
  fold_dat <- tibble(id = 1:nrow(plot_data), folds)
  fold_index <- split(fold_dat$id, fold_dat$folds)

  # Fit SuperLearner model on full data
  sl_dr <- SuperLearner(Y = y_dr_sorted,
                        X = X_train,
                        SL.library = SL.library,
                        family = gaussian(),
                        cvControl = list(V = num.folds, validRows = fold_index))

  # Create prediction grid
  x_grid <- seq(min(x_sorted), max(x_sorted), length.out = 200)
  X_pred <- data.frame(x = x_grid)

  # Get predictions
  pred_dr <- predict(sl_dr, newdata = X_pred, onlySL = TRUE)$pred

  # Create prediction dataframe
  pred_data <- data.frame(
    log_var = x_grid,
    pred_dr = pred_dr
  )

  # Keep original data for rug plot
  plot_data$log_var <- x_vals

  # Create plot
  p <- ggplot() +
    geom_line(data = pred_data, aes(x = log_var, y = pred_dr),
              color = "black", linewidth = 1) +
    geom_rug(data = plot_data, aes(x = log_var), sides = "b", length = unit(2, "mm")) +
    geom_hline(yintercept = ate_estimates, color = "gray", linetype = "dashed") +
    xlab(x_label) + ylim(-8, -2) +
    theme_classic(base_size = 20)

  # Conditionally add y-axis label
  if (show_ylab) {
    p <- p + ylab("Age-AFC Association")
  } else {
    p <- p + ylab(NULL)
  }

  # Print SuperLearner weights for diagnostic purposes
  cat("\nSuperLearner weights for", var_name, "(DR scores):\n")
  print(sl_dr$coef)
  cat("\n")

  return(p)
}

# Create plots for all 17 EDCs using SuperLearner
var_names <- c("MCOP","MiBP", "MEHP")
ncol_grid <- 3
plot_list <- lapply(seq_along(var_names), function(i) {
  var <- var_names[i]
  show_ylab <- (i - 1) %% ncol_grid == 0
  create_cate_plot_sl(var, afc_clean_notrunc, x_labels[var], res_dat$estimate,
                      show_ylab = show_ylab)
})

# Combine plots
plot_grid <- grid.arrange(grobs = plot_list, ncol = ncol_grid)

ggsave(plot = plot_grid, filename = here("figures", "cate_functions_paper.png"),
       width = 20, height = 16, units = "cm", dpi = 300)
