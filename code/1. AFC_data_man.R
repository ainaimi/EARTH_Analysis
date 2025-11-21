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
load(here("data", "imputed_EARTH.Rdata")) 
new_afc <- afc_clean <- imputed_data %>% 
  rename(AFCt = AFC_t)

## some data basics, names, dimensions, etc
names(afc_clean) 

## A question for Audrey: so what is DOR?
afc_summary <- afc_clean %>% 
  group_by(DOR) %>% 
  summarize(n(),
            min(AFCt),
            max(AFCt),
            mean(AFCt),
            median(AFCt))

saveRDS(afc_summary, here("output", "afc_summary_by_DOR.rds"))

# the relevant data:
the_data <- c("AFCt", "year", "month", "DOR","eversmok",
              "age", "bmi", "races", "educ1", "smokstat", "previousIVF", "previousIUI", "gravid",

              # 17 EDC variables (16 SG-adjusted + 1 Hg)
              "MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "sumDEHP",
              "BPA", "BP", "MP", "PP",
              "Hg")

# Non-redundant imputation flags (linearly dependent flags removed)
impute_flags <- c("imp_sgratio_pht", "imp_smokstat", "imp_races", "imp_bmi", "imp_previousIUI",
                  "imp_age", "imp_educ1", "imp_Hg", "imp_mBP", "imp_mCNP",
                  "imp_B_PB", "imp_AFScanDate")

# 17 EDC variables with <40% missing (16 SG-adjusted + 1 Hg)
env_vars <- c("MBP", "MiBP", "MCNP", "MCOP", "MECPP", "MEHHP", "MEHP", "MEOHP",
              "MCPP", "MEP", "MBzP", "sumDEHP",
              "BPA", "BP", "MP", "PP",
              "Hg")

length(env_vars)

skim_results <- skimr::skim(new_afc[,the_data])

saveRDS(skim_results, "skim_results.rds")

summary_stats <- new_afc %>%
  summarise(across(c(AFCt, year, month, age,
                     MBP, MiBP, MCNP, MCOP, MECPP, MEHHP, MEHP, MEOHP,
                     MCPP, MEP, MBzP, sumDEHP,
                     BPA, BP, MP, PP, Hg), .fns = 
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

afc_clean_notrunc <- afc_clean[,c(the_data, impute_flags)]
afc_clean_trunc   <- afc_clean[,c(the_data, impute_flags)]
afc_clean_trunc[,env_vars] <- apply(afc_clean_trunc[,env_vars], 2, trunc_func)

skim(afc_clean)
skim(afc_clean_trunc)

save(afc_clean_trunc, 
     file = here("data", "afc_clean_trunc.Rdata"))

save(afc_clean_notrunc, 
     file = here("data", "afc_clean_notrunc.Rdata"))