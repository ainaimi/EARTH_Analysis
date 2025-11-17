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
  mgcv
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

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg",#"miBP_sg",
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
  geom_point(aes(x = abs(statistic.x), 
                 y = abs(statistic.y))) +
  ylim(0,2) + xlim(0,2) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
  geom_text(aes(label = term, x = abs(statistic.x), y = abs(statistic.y)), vjust = -0.5) +
  xlab("Test Statistic, DR Learner") + 
  ylab("Test Statistic, Causal Forest") +
  theme_classic(base_size = 16)

ggsave(filename = here("figures", "teststat_scatter.pdf"), 
       width = 16, height = 16, units = "cm")


## select top four environmental variables
var_names <- table2[order(-abs(table2$statistic.x)), ]$term[1:4]

env_if <- data.frame(
  melt(afc_clean_notrunc[,var_names]),
  forest_CATE = rep(if_data$forest_scores, length(var_names)),
  dr_learner_CATE = rep(if_data$dr_scores, length(var_names))
)

ggplot(env_if) +
  geom_smooth(aes(x = log(value), y = forest_CATE),
              method = "lm", se = F, color = "black") +
  geom_smooth(aes(x = log(value), y = dr_learner_CATE),
              method = "lm", se = F, color = "black", linetype = "dashed") +
  geom_rug(aes(x = jitter(log(value), 3)),
           length = unit(2, "mm")) +
  geom_hline(yintercept = res_dat$estimate, color = "gray") +
  facet_wrap(~variable, scale = "free") + ylab("Age-AFC Association") + xlab("log(Value)") +
  theme_classic(base_size = 16)

ggsave(filename = here("figures", "cate_functions.pdf"), 
       width = 20, height = 20, units = "cm")