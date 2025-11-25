pacman::p_load(
  rio,
  here,
  tidyverse,
  parallel,
  SuperLearner,
  broom,
  vip,
  ranger,
  xgboost,
  polspline,
  earth
)

load(here("data", "afc_clean_notrunc.Rdata"))

afc_clean_notrunc$age_bin <- as.numeric(afc_clean_notrunc[,"age"] >= 35)

ggplot(afc_clean_notrunc) +
  geom_point(aes(x = age, y = AFCt)) +
  geom_smooth(aes(x = age, y = AFCt))

ggsave(filename = here("figures", "main_scatterplot.png"), dpi = 300)

descriptives_table1 <- afc_clean_notrunc %>% 
  group_by(as.numeric(age>=35)) %>% 
  summarize(mean_bmi = mean(bmi),
            mean_gravid = mean(as.numeric(gravid)-1),
            mean_IUI = mean(as.numeric(previousIUI)-1),
            mean_IVF = mean(as.numeric(previousIVF)-1),
            mean_AFC = mean(AFCt))

saveRDS(descriptives_table1, here("output", "descriptives_table1.rds"))

# create exposure variable
exposure <- as.numeric(unlist(afc_clean_notrunc[,"age_bin"]))

table(exposure)
mean(exposure)

# create outcome variable
outcome <- as.numeric(unlist(afc_clean_notrunc[,"AFCt"]))
hist(outcome)

# identify and transform categorical
categorical_vars <- afc_clean_notrunc %>% select(where(is.factor), -DOR)

# 17 EDC variables with <40% missing (16 SG-adjusted + 1 Hg)
# (defined in file 1. AFC_data_man.R)
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", #"sumDEHP",
              "BPA", "BP", "MP", "PP",
              "Hg")

# identify and transform continuous
continuous_vars <- afc_clean_notrunc %>%
  select(where(is.numeric), -starts_with("imp_"), -age_bin, -AFCt, -age, -DOR)

# identify imputation flags
flag_vars <- afc_clean_notrunc %>% select(starts_with("imp_"))

new_afc <- data.frame(outcome, exposure, categorical_vars, continuous_vars, flag_vars)

# relevant covariate set
covariates <- c(names(continuous_vars), names(categorical_vars), names(flag_vars))

# use covariate set to construct design matrix that doesn't include exposure
covariates_matrix <- data.frame(model.matrix(formula(paste0("~", 
                                                            paste0(covariates, collapse="+"))), 
                                             data=new_afc)[,-1])

# use covariate set to construct design matrix that includes exposure
covariates_matrix_w <- data.frame(covariates_matrix, age_bin = exposure)

# for reproducibility, we set the seed
set.seed(123)

# construct a sample size variable
n <- nrow(new_afc)

## DR learner with SuperLearner
mean_learner <- "SL.mean"
glm_learner <- "SL.glm"

# ranger learner
ranger_learner <- create.Learner("SL.ranger",
                                 params = list(min.node.size = 25),
                                 tune = list(num.trees = c(500,1000),
                                             mtry = c(3,4,5)))

# glmnet learner
glmnet_learner <- create.Learner("SL.glmnet",
                                 tune = list(alpha = seq(0,1,.2)))



# xgboost learner
xgboost_learner <- create.Learner(base_learner = "SL.xgboost",
                                  params = list(minobspernode = 25),
                                  tune = list(
                                    max_depth = c(2, 4),
                                    nrounds = c(500, 1000)
                                            )
                                  )

sl_lib <- list("SL.mean", "SL.glm", 
               "SL.ranger_1", "SL.ranger_2", "SL.ranger_3", "SL.ranger_4", "SL.ranger_5", "SL.ranger_6",
               "SL.xgboost_1", "SL.xgboost_2", "SL.xgboost_3", "SL.xgboost_4",
               "SL.glmnet_1", "SL.glmnet_2", "SL.glmnet_3", "SL.glmnet_4", "SL.glmnet_5", "SL.glmnet_6",
               c("SL.glm.interaction", "screen.corP"))
  

# Specify the number of folds for V-fold cross-validation
# Use same folds as used for causal_forest function
# Doing cross-validation this way automatically deploys cross-fitting
num.folds <- 10
folds <- sort(seq(length(exposure)) %% num.folds) + 1
fold_dat <- tibble(id = 1:length(exposure),folds)
fold_index <- split(fold_dat$id,fold_dat$folds)

# we'll use parallel processing to speed things up
options(mc.cores = detectCores() - 2)

# covariates_matrix_w_augment <- cbind(covariates_matrix_w, augment_data)

if(file.exists(here("misc","fit_mu.RDS"))){
  fit_mu <- readRDS(here("misc","fit_mu.RDS"))
} else{
  fit_mu <- CV.SuperLearner(Y = outcome,
                            X = covariates_matrix_w, 
                            method = "method.NNLS", 
                            family = gaussian(),
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),
                            control = list(saveCVFitLibrary = T),
                            parallel = "multicore",
                            verbose = T)
  saveRDS(object = fit_mu, file = here("misc","fit_mu.RDS"))
}

if(file.exists(here("misc","fit_pi.RDS"))){
  fit_pi <- readRDS(here("misc","fit_pi.RDS"))
} else{
  fit_pi <- CV.SuperLearner(Y = exposure,
                            X = covariates_matrix,
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),#, stratifyCV = TRUE),
                            control = list(saveCVFitLibrary = F),
                            parallel = "multicore",
                            verbose = T)
  saveRDS(object = fit_pi, file = here("misc","fit_pi.RDS"))
}

summary(fit_mu)

saveRDS(summary(fit_mu), here("output", "mu_sl_summary.rds"))

summary(fit_pi)

saveRDS(summary(fit_pi), here("output", "pi_sl_summary.rds"))

## generate variable importance
pred_wrapper <- function(sl, newdata) {
  as.numeric(predict(sl, newdata, onlySL=T)$pred)
}

fit_mu_2 <- SuperLearner(Y = outcome,
                         X = covariates_matrix,
                         method = "method.NNLS", 
                         family = gaussian(),
                         SL.library = list("SL.mean", "SL.glm", 
                                           "SL.ranger_1", "SL.ranger_2", "SL.ranger_3", "SL.ranger_4", "SL.ranger_5", "SL.ranger_6",
                                           "SL.xgboost_1", "SL.xgboost_2", "SL.xgboost_3", "SL.xgboost_4",
                                           "SL.glmnet_1", "SL.glmnet_2", "SL.glmnet_3", "SL.glmnet_4", "SL.glmnet_5", "SL.glmnet_6",
                                           c("SL.glm.interaction", "screen.corP")),
                         cvControl = list(V = num.folds, validRows = fold_index),
                         control = list(saveCVFitLibrary = T),
                         verbose = T)

VarImpSL <- vi_permute(object = fit_mu_2, 
                       #method = "permute", 
                       feature_names = names(covariates_matrix), 
                       train = covariates_matrix, 
                       target = as.numeric(outcome),
                       metric = "RMSE",
                       smaller_is_better = T,
                       type = "difference", 
                       nsim = 5, 
                       pred_wrapper = pred_wrapper) 

VarImpSL <- VarImpSL[,c("Variable","Importance")]

VarImpSL$Model <- "Super Learner"

VarImp <- VarImpSL

VarImp <- VarImp %>% arrange(desc(Importance))
var_names <- unique(VarImp$Variable)

ggplot(VarImp) +
  geom_bar(aes(x = reorder(Variable, -Importance),
               y = Importance),
           fill = "steelblue",
           stat="identity") +
  theme(axis.text.x = element_text(angle = 70,
                                   vjust = 0.5,
                                   hjust=.5)) +
  scale_y_continuous(expand = c(0,0))

ggsave(filename = here("figures", "variable_importance_plot.png"), units = "cm", height = 16, width = 16, dpi = 300)

## create variable list for environmental chemicals
var_names_sl <- unique(VarImpSL[order(-VarImpSL$Importance),]$Variable)
var_names_sl <- var_names_sl[var_names_sl %in% env_vars] # to save for later...

# ps overlap plot
plot_dat_ps <- tibble(`Propensity Score` = fit_pi$SL.predict,
                      Exposure = factor(exposure))

ps_overlap_plot <- ggplot(plot_dat_ps) +
  geom_density(aes(x = `Propensity Score`,
                   group = Exposure,
                   fill = Exposure), alpha = .5) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Propensity Score Overlap")

ps_overlap_plot

ggsave(plot = ps_overlap_plot, filename = here("figures", "ps_overlap_plot.png"),
       units = "cm", height = 16, width = 16, dpi = 300)

# Again, with these diagnostics complete, we can construct and ATE estimate using the double robust AIPW estimator, with the exposure and outcome models obtained via stacking:

## cross-fit predictions

pscore <- as.matrix(fit_pi$SL.predict)

mu_hat <- as.matrix(fit_mu$SL.predict)

mu_hat1 <- NULL
for(i in 1:num.folds){
  mu_hat1 <- rbind(mu_hat1, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = base::transform(
                             covariates_matrix_w[fold_index[[i]],], age_bin = 1), 
                           onlySL=T)$pred)
}

mu_hat0 <- NULL
for(i in 1:num.folds){
  mu_hat0 <- rbind(mu_hat0, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = base::transform(
                             covariates_matrix_w[fold_index[[i]],], age_bin = 0), 
                           onlySL=T)$pred)
}

## aipw
aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1){
  aipw_score <- ((2*exposure - 1)*(outcome - mu_hat))/((2*exposure - 1)*pscore + (1 - exposure)) + (mu_hat1 - mu_hat0)
  return(aipw_score)
}

aipw_score <- aipw_func(exposure, 
            outcome, 
            pscore, 
            mu_hat, 
            mu_hat0, 
            mu_hat1)


names(aipw_score) <- NULL

aipw_psi <- mean(aipw_score)

aipw_se <- sd(aipw_score)/sqrt(n)

aipw_ate <- c(aipw_psi, aipw_se)

aipw_ate <- cbind(t(aipw_ate), 
                  aipw_ate[1] - 1.96*aipw_ate[2],
                  aipw_ate[1] + 1.96*aipw_ate[2])

colnames(aipw_ate) <- c("estimate", "std.err", "LCL", "UCL")

aipw_ate

res_ate <- aipw_ate

row.names(res_ate) <- c("AIPW with Stacking")

res_ate

saveRDS(res_ate, here("output", "ate_results.rds"))

afc_clean_notrunc <- afc_clean_notrunc %>%
  mutate(dr_scores = unlist(aipw_score))

save(afc_clean_notrunc, 
     file = here("data", "afc_clean_notrunc_IF.Rdata"))