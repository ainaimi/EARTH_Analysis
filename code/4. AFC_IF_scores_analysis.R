pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  reshape2,
  clubSandwich
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

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

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


lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>% 
  select(-forest_scores) %>% 
  mutate(across(all_of(env_vars), log))
dr_learner <- lm(dr_scores ~ ., data = lin_proj)
tidy(dr_learner)

lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>% 
  select(-dr_scores) %>% 
  mutate(across(all_of(env_vars), log))
causal_for <- lm(forest_scores ~ ., data = lin_proj)
tidy(causal_for)

table2 <- left_join(tidy(dr_learner), 
          tidy(causal_for), 
          by = "term") %>% 
  filter(term != "(Intercept)") %>% 
  select(term, estimate.x, std.error.x, statistic.x, estimate.y, std.error.y, statistic.y)

ggplot(table2) +
  geom_point(aes(x = statistic.x, 
                 y = statistic.y)) +
  ylim(-2,2) + xlim(-2,2) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
  geom_text(aes(label = term, x = statistic.x, y = statistic.y), vjust = -0.5) +
  xlab("Test Statistic, DR Learner") + 
  ylab("Test Statistic, Causal Forest") +
  theme_classic(base_size = 16)

ggsave(filename = here("figures", "teststat_scatter.pdf"), 
       width = 16, height = 16, units = "cm")


library(mgcViz)
library(mgcv)
lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>% 
  select(-forest_scores)  %>% 
  mutate(across(all_of(env_vars), log))
dr_learner_gam <- gam(dr_scores ~ s(dehp_sg, bs = "tp", k = 15) ,
                   method = "REML", data = lin_proj)
summary(dr_learner_gam)
b <- getViz(dr_learner_gam)
# print() only needed because we want to plot on a single page
print(plot(b), pages = 1)

lin_proj <- data.frame(if_data, afc_clean_notrunc[,env_vars]) %>% 
  select(-dr_scores)  %>% 
  mutate(across(all_of(env_vars), log))
dr_learner_gam <- gam(forest_scores ~ s(mEP1_sg) + s(mBZP1_sg) + s(mCNP_sg) + s(miBP_sg) +
                        s(mBP_sg) + s(BP_3_sg) + s(M_PB_sg) + s(dehp_sg) + s(BPA_sg) +
                        s(mCOP_sg) + s(mCPP_sg) + s(B_PB_sg) + s(P_PB_sg) + s(Hg), data = lin_proj)
summary(dr_learner_gam)
plot(dr_learner_gam)

# if plots
ggplot(if_data) +
  geom_point(aes(x = forest_scores, 
                 y = dr_scores))

ifM <- melt(if_data)

ggplot(ifM, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable)

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

env_if <- data.frame(
  melt(afc_clean_notrunc[,env_vars]),
  forest_CATE = rep(if_data$forest_scores, length(env_vars)),
  dr_learner_CATE = rep(if_data$dr_scores, length(env_vars))
)

# merge_vars <- env_if %>% 
#   filter(variable %in% c("BPA_sg", "Hg")) %>% 
#   group_by(variable) %>%
#   summarize(q1 = quantile(value, 0.1),
#             q3 = quantile(value, 0.9))
# 
# merge_vars
# 
# env_if <- merge(env_if, merge_vars, by = "variable")

head(env_if)

env_if %>% #[env_if$value>=env_if$q1&env_if$value<=env_if$q3,] %>% 
  ggplot(.) +
  # geom_point(aes(x = log(value), y = forest_CATE),
  #             color = "gray") +
  geom_smooth(aes(x = log(value), y = forest_CATE),
              method = "lm", se = F, color = "black") +
  geom_smooth(aes(x = log(value), y = dr_learner_CATE),
              method = "lm", se = F, color = "black", linetype = "dashed") +
  geom_smooth(aes(x = log(value), y = forest_CATE),
              method = "loess", span = 1, se = F, color = "gray") +
  geom_smooth(aes(x = log(value), y = dr_learner_CATE),
              method = "loess", span = 1, se = F, color = "gray", linetype = "dashed") +
  geom_rug(aes(x = jitter(log(value), 3)),
           length = unit(2, "mm")) +
  # ylim(-.5,.5) +
  geom_hline(yintercept = res_dat$estimate, color = "red") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  facet_wrap(~variable)


afc_clean_notrunc %>% 
  ggplot(.) + geom_point(aes(x = jitter(age, factor = 2), y = AFC_t ))


## lepski for dr learner
# 
# n <- nrow(env_if[env_if$variable == "BP_3_sg",])
# X <- env_if[env_if$variable == "BP_3_sg",]$value
# Y <- env_if[env_if$variable == "BP_3_sg",]$dr_learner_CATE
# 
# res <- data.frame(X,Y)
# 
# res <- res[order(res$X),]
# 
# X <- log(res$X)
# Y <- res$Y
# 
# # Define a range of bandwidths
# bandwidths <- seq(0.01, 2, length.out = 60)
# n_bandwidths <- length(bandwidths)
# 
# lepski_kernel_regression <- function(X, Y, bandwidths, kernel = "gaussian", lambda = 1.5) {
#   
#   n <- length(X)
#   n_bandwidths <- length(bandwidths)
#   f_hat <- matrix(NA, nrow = n, ncol = n_bandwidths)
#   
#   # Define kernel function
#   kernel_function <- function(u) {
#     if (kernel == "gaussian") {
#       return(dnorm(u))
#     } else if (kernel == "epanechnikov") {
#       return(0.75 * (1 - u^2) * (abs(u) <= 1))
#     } else {
#       stop("Unsupported kernel. Choose 'gaussian' or 'epanechnikov'.")
#     }
#   }
#   
#   # Kernel Regression
#   for (j in 1:n_bandwidths) {
#     h <- bandwidths[j]
#     for (i in 1:n) {
#       weights <- kernel_function((X - X[i]) / h)
#       f_hat[i, j] <- sum(weights * Y) / sum(weights)
#     }
#   }
#   
#   # Lepski’s method to select the best bandwidth
#   lepski_method <- function(f_hat, lambda) {
#     discrepancies <- matrix(NA, nrow = n, ncol = n_bandwidths)
#     
#     for (j in 1:n_bandwidths) {
#       for (k in j:n_bandwidths) {
#         if (j == k) {
#           discrepancies[, j] <- 0
#         } else {
#           discrepancies[, j] <- apply(f_hat[, j:k], 1, function(row) {
#             max(abs(row - row[1]))
#           })
#         }
#       }
#     }
#     
#     selected_bandwidth_index <- rep(NA, n)
#     for (i in 1:n) {
#       for (j in 1:n_bandwidths) {
#         if (all(discrepancies[i, j:n_bandwidths] <= lambda * sd(Y))) {
#           selected_bandwidth_index[i] <- j
#           break
#         }
#       }
#     }
#     
#     return(selected_bandwidth_index)
#   }
#   
#   selected_indices <- lepski_method(f_hat, lambda)
#   selected_bandwidths <- bandwidths[na.omit(selected_indices)]
#   optimal_index <- median(na.omit(selected_indices))
#   
#   return(list(
#     estimate = f_hat[, optimal_index],
#     bandwidth = bandwidths[optimal_index],
#     f_hat = f_hat,
#     selected_indices = selected_indices
#   ))
# }
# 
# 
# result <- lepski_kernel_regression(X, Y, bandwidths, kernel = "gaussian", lambda = 1.5)
# 
# # Plot results
# plot(X, Y, pch = 16, col = "blue", main = "Lepski's Method - Selected Bandwidth")
# lines(X, result$estimate, col = "darkgreen", lwd = 2)
# 
# 

# fit your model (replace df / outcome / predictors accordingly)
fit <- lm(outcome ~ mEP1_sg + mBZP1_sg + mCNP_sg + miBP_sg + mBP_sg + BP_3_sg +
            M_PB_sg + dehp_sg + BPA_sg + mCOP_sg + mCPP_sg + B_PB_sg +
            P_PB_sg + Hg,
          data = df)

# ---- Per-term small-sample robust results (CR2) ----
# coef_test gives small-sample adjusted t-tests; we convert to chi-square (df=1)
ct <- clubSandwich::coef_test(
  dr_learner,
  vcov = "CR2",           # small-sample adjusted sandwich
  test = "Satterthwaite"  # small-sample df for t-tests
)

term_tbl <- ct %>%
  transmute(
    term        = term,
    estimate    = beta,
    std.error   = SE,
    statistic   = tstat,                             # robust t
    p.value     = p_Satt,                            # Satterthwaite p
    chisq_rob   = tstat^2,                           # Wald chi-square (df = 1)
    p_chisq_rob = pchisq(chisq_rob, df = 1, lower.tail = FALSE)
  )

# ---- Global robust Wald χ²: H0: all slopes = 0 ----
# Build constraint matrix selecting all non-intercept coefficients
B <- diag(length(coef(fit)))[-1, , drop = FALSE]         # drop intercept row
colnames(B) <- names(coef(fit))

global_test <- clubSandwich::Wald_test(
  fit,
  constraints = B,
  vcov = "CR2",
  test = "chi2"   # returns a chi-square with small-sample scaling
)

# View results
term_tbl
global_test
