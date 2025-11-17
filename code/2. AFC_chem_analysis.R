pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  reshape2
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

## goal of this program is to explore distribution of environmental exposures 
## in EARTH data
## 
## Compare both trimmed and untrimmed versions

load(here("data", "afc_clean_trunc.Rdata")) 

load(here("data", "afc_clean_notrunc.Rdata")) 

head(afc_clean_notrunc)
head(afc_clean_trunc)

## some basic outlier / leverage analyses

## do the people with extreme chem values also have large outcome values?

env_vars <- c("mEP1_sg", "mBZP1_sg", "mCNP_sg","miBP_sg",
              "mBP_sg","BP_3_sg", "M_PB_sg", "dehp_sg", 
              "BPA_sg", "mCOP_sg", "mCPP_sg", "B_PB_sg", "P_PB_sg", "Hg")

dim(afc_clean_trunc[,env_vars])

formulaVars <- paste(env_vars, collapse = "+")
modelForm <- as.formula(paste0("AFC_t ~", formulaVars))
modelForm

mod <- lm(modelForm, data = afc_clean_notrunc)

plot(mod)

diagnostic_data <- broom::augment(mod)

# problem_row <- c(116, 117, 194, 195, 289, 290, 620, 621)
# 
# diagnostic_data[problem_row,]

dd <- diagnostic_data[,16:21]

ddM <- melt(dd)

ggplot(ddM, aes(x=value)) +
  geom_histogram() + 
  facet_wrap(~variable, scales = "free")



