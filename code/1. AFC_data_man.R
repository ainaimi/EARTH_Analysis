pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse
)

#Load afc data, double check variables + attributes
load(here("data", "imputed_EARTH.Rdata")) 
afc_clean <- imputed_data %>%
  rename(AFCt = AFC_t) 

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

skim_results <- skimr::skim(afc_clean[,the_data])

saveRDS(skim_results, here("output", "skim_results.rds"))

summary_stats <- afc_clean %>%
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

saveRDS(summary_stats, here("output", "summary_stats.rds"))

# trim the data
trunc_func <- function(x, upper_bound = .925){
  x[x > quantile(x, upper_bound)] <- quantile(x, upper_bound)
  return(x)
}

afc_clean_notrunc <- afc_clean[,c(the_data, impute_flags)]
afc_clean_trunc   <- afc_clean[,c(the_data, impute_flags)]
afc_clean_trunc[,env_vars] <- apply(afc_clean_trunc[,env_vars], 2, trunc_func)

save(afc_clean_trunc, 
     file = here("data", "afc_clean_trunc.Rdata"))

save(afc_clean_notrunc, 
     file = here("data", "afc_clean_notrunc.Rdata"))