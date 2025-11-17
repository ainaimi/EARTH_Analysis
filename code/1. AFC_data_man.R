pacman::p_load(
  rio,          
  here,         
  skimr,    
  tidyverse,
  parallel,
  scales,
  haven,
  lmtest,
  sandwich,
  grf,
  SuperLearner,
  pROC,
  broom,
  xtable,
  gridExtra, 
  vip,
  ranger,
  xgboost,
  polspline,
  earth
)

# set the theme for figures
thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    text = element_text(size = 20)
  )
theme_set(thm)

#Load afc data, double check variables + attributes
afc_clean <- read_csv(here("data","afc_clean.csv"))
new_afc <- read_csv(here("data","afc_clean.csv"))

# afc_clean <- afc_clean %>% mutate(age_bin = as.numeric(age>=35))
# new_afc <- new_afc %>% mutate(age_bin = as.numeric(age>=35))

names(afc_clean)

names(new_afc)

## so what is DOR?
afc_clean %>% 
  group_by(DOR) %>% 
  summarize(n(),
            min(AFC_t),
            max(AFC_t),
            mean(AFC_t),
            median(AFC_t))

# the relevant data:
the_data <- c("AFC_t", "year", "DOR",
              "age", "bmi", "white", "educ1", "eversmok", "previousIVF", "previousIUI", "gravid",
              "mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

skimr::skim(new_afc[,the_data])

afc <- afc_clean %>% rename(PPB     = P_PB, 
                            PPBsg   = P_PB_sg,
                            mCNPsg  = mCNP_sg,
                            mBPsg   = mBP_sg,
                            AFCt    = AFC_t,
                            mEP1sg  = mEP1_sg, 
                            mBZP1sg = mBZP1_sg,
                            miBPsg  = miBP_sg,
                            BP3sg   = BP_3_sg, 
                            MPBsg   = M_PB_sg, 
                            dehpsg  = dehp_sg,
                            BPAsg   = BPA_sg, 
                            mCOPsg  = mCOP_sg, 
                            mCPPsg  = mCPP_sg, 
                            BPBsg   = B_PB_sg)

afc %>% summarise(across(c(AFCt, year, PPB, PPBsg, mCNP, mCNPsg, mBP, mBPsg, Hg,
                           PPB, PPBsg, mCNPsg, mBPsg, AFCt, mEP1sg, mBZP1sg, miBPsg,
                           BP3sg, MPBsg, dehpsg, BPAsg, mCOPsg, mCPPsg, BPBsg), .fns = 
                           list(min = min,
                                median = median,
                                mean = mean,
                                stdev = sd,
                                q25 = ~quantile(., 0.25),
                                q75 = ~quantile(., 0.75),
                                max = max))) %>% 
  round(., digits = 2) %>% 
  pivot_longer(everything(), names_sep="_", names_to=c('variable', '.value'))

trunc_func <- function(x, upper_bound = .925){
  x[x > quantile(x, upper_bound)] <- quantile(x, upper_bound)
  return(x)
}

afc_clean_notrunc <- afc_clean[,the_data]
afc_clean_trunc <- afc_clean[,the_data]
afc_clean_trunc[,env_vars] <- apply(afc_clean_trunc[,env_vars], 2, trunc_func)

skim(afc_clean)
skim(afc_clean_trunc)



save(afc_clean_trunc, 
     file = here("data", "afc_clean_trunc.Rdata"))

save(afc_clean_notrunc, 
     file = here("data", "afc_clean_notrunc.Rdata"))


##### STOP


# Test to see if the ATEs from causal forest versus AIPW are different:

diff_eff <- aipw_ate[1] - cf_ate[1]

se_diff <- sqrt(var(aipw_score) + var(forest_score) - 2*cov(aipw_score, forest_score))/sqrt(n)

diff_eff

lcl <- diff_eff - 1.96*se_diff
ucl <- diff_eff + 1.96*se_diff

c(diff_eff, lcl, ucl)

2*pnorm(abs(diff_eff/se_diff), lower.tail=FALSE)
2*pt( abs(diff_eff/se_diff), df = n-1, lower.tail=FALSE)

## CATEs with causal forest

tau.hat = predict(forest)$predictions

# Compare regions with high and low estimated CATEs
high_effect = tau.hat > median (tau.hat)
ate.high = average_treatment_effect(forest, subset = high_effect)
ate.low = average_treatment_effect(forest, subset = !high_effect)

ate.high
ate.low

# Run best linear predictor analysis
test_calibration(forest)

blp_forest <- best_linear_projection(forest, covariates_matrix[,c(var_names_rf)])

score_dat <- data.frame(aipw_score, covariates_matrix[,c(var_names_sl)])
score_fit <- lm(aipw_score ~ ., data = score_dat)

blp_aipw <- coeftest(score_fit, vcov = vcovHC(score_fit, type = "HC3"))

blp_forest <- tidy(blp_forest)[,c(1,2,3,4,5)]
blp_aipw <- tidy(blp_aipw)[,c(1,2,3,4,5)]

names(blp_forest) <- c("term", "est.forest", "se.forest", "test.stat.forest", "p.forest")
names(blp_aipw) <- c("term", "est.aipw", "se.aipw", "test.stat.aipw", "p.aipw")

blp_res <- left_join(blp_forest,
                     blp_aipw,
                     by = "term")

blp_res <- blp_res %>% 
  filter(term != "(Intercept)") %>% 
  arrange(desc(abs((test.stat.forest+test.stat.aipw)/2)))


## super learner based modifier plots
## 

## change library to more regression based

sl_lib <- list("SL.mean","SL.earth","SL.gam","SL.glm","SL.loess","SL.polymars")

num.folds <- 50
folds <- sort(seq(n) %% num.folds) + 1
fold_dat <- tibble(id = 1:n,folds)
fold_index <- split(fold_dat$id,fold_dat$folds)

y_CF <- (covariates_matrix$forest_score)
y_DR <- (covariates_matrix$aipw_score)
x <- covariates_matrix %>% select(mCPP_sg)

CF <- CV.SuperLearner(Y = y_CF,
                      X = x, 
                      method = "method.NNLS", 
                      family = gaussian,
                      SL.library = sl_lib,
                      cvControl = list(V = num.folds, validRows = fold_index),
                      control = list(saveCVFitLibrary = T),
                      parallel = "seq",
                      verbose = T)

CF$coef

summary(CF)

DR <- CV.SuperLearner(Y = y_DR,
                      X = x, 
                      method = "method.NNLS", 
                      family = gaussian,
                      SL.library = sl_lib,
                      cvControl = list(V = num.folds, validRows = fold_index),
                      control = list(saveCVFitLibrary = T),
                      parallel = "seq",
                      verbose = T)

#Generate the predicted values 
m.value <- data.frame(seq(min(x), max(x), by = .01))
colnames(m.value)[1] <- c("mCPP_sg")

cf.value<- NULL
for(i in 1:num.folds){
  cf.value <- rbind(cf.value, 
                    predict(CF$AllSL[[i]],
                            newdata = m.value, 
                             onlySL=T)$pred)
}
cf.value_f<-data.frame(cbind(m.value,cf.value))
colnames(cf.value_f)[2] <- c("mean_difference")
cf.value_f<- cf.value_f %>%
  group_by(mCPP_sg) %>%
  summarize(mean_difference = mean(mean_difference))
#cf.value_f$Hg <- cf.value_f$Hg*(max(afc_clean$Hg) - min(afc_clean$Hg)) + min(afc_clean$Hg)
cf.value_f$method<-c("Causal Forest")

dr.value<- NULL
for(i in 1:num.folds){
  dr.value <- rbind(dr.value, 
                    predict(DR$AllSL[[i]],
                            newdata = m.value, 
                            onlySL=T)$pred)
}
dr.value_f<-data.frame(cbind(m.value,dr.value))
colnames(dr.value_f)[2] <- c("mean_difference")
dr.value_f<- dr.value_f %>%
  group_by(mCPP_sg) %>%
  summarize(mean_difference = mean(mean_difference))
#dr.value_f$Hg <- dr.value_f$Hg*(max(afc_clean$Hg) - min(afc_clean$Hg)) + min(afc_clean$Hg)
dr.value_f$method<-c("DR Learner")

##Plot
plot_dat <- rbind(cf.value_f,dr.value_f)
write.csv(plot_dat, here("output","Data_Hg_Plot_main.csv"))

score1 <- score_dat %>% select(aipw_score, mCPP_sg)

score1$method <- "DR Learner"

score2 <- data.frame(score = forest_score, mCPP_sg = afc_clean_trunc$mCPP_sg, method = "Causal Forest")

names(score1)[1] <- "score"

score_dat_trans <- rbind(score1, score2)

f1 <- ggplot() +
  geom_rug(data = afc_clean_trunc, 
           aes(x = jitter(mCPP_sg, 3)),
           length = unit(2, "mm")) +
  geom_point(data = score_dat_trans,
            aes(y = score,
                x = mCPP_sg,
                group = method,
                color = method)) +
  geom_line(data = plot_dat,
            aes(y = mean_difference,
                x = mCPP_sg,
                group = method,
                color = method)) +
  scale_y_continuous(expand = c(.01,.01)) +
  scale_x_continuous(expand = c(.01,.01)) +
  xlab("mCPP_sg") +
  ylab("Mean Difference")

f1



