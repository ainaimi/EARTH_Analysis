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

## some data basics, names, dimensions, etc
names(afc_clean) 

## A question for Audrey: so what is DOR?
afc_summary <- afc_clean %>% 
  group_by(DOR) %>% 
  summarize(n(),
            min(AFC_t),
            max(AFC_t),
            mean(AFC_t),
            median(AFC_t))

saveRDS(afc_summary, here("output", "afc_summary_by_DOR.rds"))

# the relevant data:
the_data <- c("AFC_t", "year", "DOR",
              "age", "bmi", "white", "educ1", "eversmok", "previousIVF", "previousIUI", "gravid",
              "mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

skim_results <- skimr::skim(new_afc[,the_data])

saveRDS(skim_results, "skim_results.rds")


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

summary_stats <- afc %>% summarise(across(c(AFCt, year, PPB, PPBsg, mCNP, mCNPsg, mBP, mBPsg, Hg,
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

saveRDS(summary_stats, "summary_stats.rds")

# trim the data
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