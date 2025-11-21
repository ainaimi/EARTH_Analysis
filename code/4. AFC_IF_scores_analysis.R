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
  ggrepel
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

load(here("data", "afc_clean_notrunc_IF.Rdata")) 

names(afc_clean_notrunc)

if_data <- afc_clean_notrunc %>% 
  select(forest_scores, dr_scores) 

if_data <- data.frame(if_data)

names(if_data) <- c("forest_scores", "dr_scores")

# 17 EDC variables with <40% missing (16 SG-adjusted + 1 Hg)
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "sumDEHP",
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

res_ate <-  rbind(ate_func(if_data[,1]),
                  ate_func(if_data[,2]))

row.names(res_ate) <- c("Causal Forest", "AIPW with Stacking")

res_ate

res_dat <- data.frame(res_ate)

# construct LaTeX Table:
# xtable(res_ate)


GGally::ggpairs(log(afc_clean_notrunc[,env_vars]))

#### DR Learner
## conditional models
lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>%
  select(-forest_scores) %>%
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


#### Causal Forest
## conditional models
lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>%
  select(-dr_scores) %>%
  mutate(across(all_of(env_vars), log))
model_cf_cond <- lm(forest_scores ~ ., data = lin_proj)
causal_for_conditional <- tidy(coeftest(model_cf_cond, vcov = vcovHC(model_cf_cond, type = "HC3"))) %>%
  filter(term != "(Intercept)")
causal_for_conditional <- causal_for_conditional  %>%
  rename_with(~ paste0("cf_", .), .cols = -term)
# unconditional models
causal_for_unconditional <- map_dfr(env_vars, function(var) {
  formula_str <- paste("forest_scores ~", var)
  model <- lm(as.formula(formula_str), data = lin_proj)
  tidy(coeftest(model, vcov = vcovHC(model, type = "HC3"))) %>%
    filter(term != "(Intercept)")  # Keep only the env_var coefficient
})
causal_for_unconditional <- causal_for_unconditional  %>%
  rename_with(~ paste0("cf_", .), .cols = -term)


table2_conditional <- left_join(dr_learner_conditional, 
                                causal_for_conditional, 
                                by = "term") %>% mutate(Type = "conditional", .before=1)


table2_unconditional <- left_join(dr_learner_unconditional, 
                                  causal_for_unconditional, 
                                  by = "term") %>% mutate(Type = "unconditional", .before=1)

table2 <- rbind(table2_unconditional, table2_conditional) 

saveRDS(table2, file = here("output", "linear_projection_output.rds"))

label_data <- table2 %>%
  filter(abs(dr_statistic) > 1.28 | abs(cf_statistic) > 1.28 | term %in% c("BP", "MBP")) %>%
  mutate(
    nudge_x_val = case_when(
      term == "Hg" & Type == "conditional" ~ 0.5,      # move right a lot
      term == "MBP" & Type == "conditional" ~ -0.5,    # move left
      term == "MBP" & Type == "unconditional" ~ -0.5,
      term %in% c("BP", "MBP", "MEHP", "MCOP") ~ 0.3,
      TRUE ~ -0.3
    ),
    nudge_y_val = case_when(
      term == "Hg" & Type == "unconditional" ~ 0.2,    # move slightly up
      term == "Hg" & Type == "conditional" ~ -0.5,     # move down a lot
      term == "BP" & Type == "conditional" ~ -0.2,     # move slightly down
      term == "MiBP" & Type == "unconditional" ~ 0.6,  # move up a lot
      term %in% c("BP", "mBP", "MEHP", "MCOP") ~ -0.3,
      TRUE ~ 0.3
    )
  )

ggplot(table2) +
geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
geom_vline(xintercept = c(-1.645, 1.645), color = "gray90", linetype = "dotted") +
geom_hline(yintercept = c(-1.645, 1.645), color = "gray90", linetype = "dotted") +
geom_vline(xintercept = c(-1.44, 1.44), color = "gray60", linetype = "dotted") +
geom_hline(yintercept = c(-1.44, 1.44), color = "gray60", linetype = "dotted") +
geom_vline(xintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
geom_hline(yintercept = c(-1.28, 1.28), color = "gray30", linetype = "dotted") +
geom_point(aes(x = dr_statistic,
               y = cf_statistic, group = Type, shape = Type), size = 3) +
  scale_shape_manual(values = c(16, 3)) +
geom_text_repel(data = label_data,
                aes(label = term, x = dr_statistic, y = cf_statistic),
                nudge_x = label_data$nudge_x_val,
                nudge_y = label_data$nudge_y_val,
                segment.color = "black",
                segment.size = 0.5,
                arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
                box.padding = 0.5,
                point.padding = 0.3) +
ylim(-3,3) + xlim(-3,3) +
xlab("Test Statistic, DR Learner") +
ylab("Test Statistic, Causal Forest") +
theme_classic(base_size = 16)+
  theme(legend.position = c(0.05, 1),
        legend.justification = c(0, 1))


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
  "sumDEHP" = "log(∑DEHP), μg/L",
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
    select(forest_scores, dr_scores, all_of(var_name)) %>%
    mutate(across(all_of(var_name), ~log(. + 0.01)))
  
  # Build formula dynamically for linear models
  formula_str <- paste0("forest_scores ~ ", var_name)
  formula_dr_str <- paste0("dr_scores ~ ", var_name)
  
  # Fit linear models
  model_forest <- lm(as.formula(formula_str), data = plot_data)
  model_dr <- lm(as.formula(formula_dr_str), data = plot_data)
  
  # Add predictions to data
  plot_data$pred_forest <- predict(model_forest)
  plot_data$pred_dr <- predict(model_dr)
  plot_data$log_var <- plot_data[[var_name]]
  
  # Create plot
  p <- ggplot(plot_data, aes(x = log_var)) +
    geom_line(aes(y = pred_forest), color = "black") +
    geom_line(aes(y = pred_dr), color = "black", linetype = "dashed") +
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
grid.arrange(grobs = plot_list, ncol = ncol_grid)

ggsave(filename = here("figures", "cate_functions_full_linear.png"),
       width = 20, height = 20, units = "cm", dpi = 300)

## subset of the EDCs
var_names <- c("MEHP", "MCOP", "MiBP")
ncol_grid <- 2
plot_list <- lapply(seq_along(var_names), function(i) {
  var <- var_names[i]
  show_ylab <- (i - 1) %% ncol_grid == 0
  create_cate_plot(var, afc_clean_notrunc, x_labels[var], res_dat$estimate,
                   show_ylab = show_ylab)
})

# Combine plots
grid.arrange(grobs = plot_list, ncol = 3)

# Function with flexible splines (GAM)
create_cate_plot <- function(var_name, data, x_label, ate_estimates, show_ylab = TRUE, k = 4) {
  # Select and clean data for this variable
  plot_data <- data %>%
    select(forest_scores, dr_scores, all_of(var_name)) %>%
    mutate(across(all_of(var_name), ~log(. + 0.01)))

  # Build formula dynamically with smoothing splines
  # k controls the basis dimension (flexibility)
  formula_str <- paste0("forest_scores ~ s(", var_name, ", k=", k, ")")
  formula_dr_str <- paste0("dr_scores ~ s(", var_name, ", k=", k, ")")

  # Fit GAM models with splines
  model_forest <- gam(as.formula(formula_str), data = plot_data)
  model_dr <- gam(as.formula(formula_dr_str), data = plot_data)

  # Add predictions to data
  plot_data$pred_forest <- predict(model_forest)
  plot_data$pred_dr <- predict(model_dr)
  plot_data$log_var <- plot_data[[var_name]]

  # Create plot
  p <- ggplot(plot_data, aes(x = log_var)) +
    geom_line(aes(y = pred_forest), color = "black") +
    geom_line(aes(y = pred_dr), color = "black", linetype = "dashed") +
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

ncol_grid <- 3
plot_list <- lapply(seq_along(var_names), function(i) {
  var <- var_names[i]
  show_ylab <- (i - 1) %% ncol_grid == 0
  create_cate_plot(var, afc_clean_notrunc, x_labels[var], res_dat$estimate,
                   show_ylab = show_ylab)
})

# Combine plots
grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = here("figures", "cate_functions_full.png"),
       width = 20, height = 20, units = "cm", dpi = 300)

## subset of the EDCs
var_names <- c("Hg", "BP")
ncol_grid <- 2
plot_list <- lapply(seq_along(var_names), function(i) {
  var <- var_names[i]
  show_ylab <- (i - 1) %% ncol_grid == 0
  create_cate_plot(var, afc_clean_notrunc, x_labels[var], res_dat$estimate,
                   show_ylab = show_ylab, k = 6)
})

# Combine plots
grid.arrange(grobs = plot_list, ncol = ncol_grid)

ggsave(filename = here("figures", "cate_functions.png"),
       width = 20, height = 20, units = "cm", dpi = 300)


























# Prepare data
plot_data <- afc_clean_notrunc %>%
  select(forest_scores, dr_scores, MEHP, MCOP, MiBP) %>%
  mutate(across(c(MEHP, MCOP, MiBP), ~log(. + 0.01))) %>%
  pivot_longer(cols = c(MEHP, MCOP, MiBP),
               names_to = "chemical",
               values_to = "log_concentration")

# Plot
ggplot(plot_data) +
  geom_point(aes(x = log_concentration, y = forest_scores),
             color = "black", alpha = 0.5) +
  geom_point(aes(x = log_concentration, y = dr_scores),
             color = "red", alpha = 0.5) +
  facet_wrap(~chemical, scales = "free_x") +
  labs(x = "Log Concentration",
       y = "IF Scores",
       title = "Forest (black) vs DR (red) Scores") +
  theme_classic()