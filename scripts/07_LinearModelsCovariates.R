#Regression of metabolites on different covariates

# Load libraries
library(tidyverse)
library(lme4)
# Source functions
source("functions/AllFunctions.R")

#Data
# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
#Add longitudinal markers: CRP and PCT and ratios and sums
data.full <- AddRatiosAndSums(AddCRPAndPCT(data.reduced))

#Covariates of interest:
#age, sex, antibiotica bij start, antibiotica switch, CURB, CRP, PCT, 
# duration of symptoms before admission, hospitalization time
#Note: "duration.of.symptoms.before.admission", "antibiotic.treatment.before.admission" have a lot of missingness

#Select variables of interest
vars_of_interest <- c("age", "sex", "curb", "hospitalization.time", "crp", "pct")

#Define metrange
metrange <- c(attr(data.full, "metrange"), attr(data.full, "ratiorange"))

#Fit a linear mixed effect model to the data for each metabolite, store estimates and t values
results_lmm <- NULL
for (i in 1:length(metrange)) {
  dat_lm <- data.full %>% 
    select(all_of(c(metrange[i], vars_of_interest, "subject.id", "day"))) %>% 
    mutate(subject.id = as.factor(subject.id))
  
  lm_model <- lmer(data = dat_lm, as.formula(paste(metrange[i], "~", paste(vars_of_interest, collapse = " + "), 
                                                   "+ day +(1 + day|subject.id)")))
  results_lmm <- rbind.data.frame(results_lmm, data.frame(metabolite = metrange[i], summary(lm_model)$coefficients[-1, ]) %>% rownames_to_column("covariate"))
}

#Fit a linear model with day30-day0 for PCT, CRP and the metabolites
#Transform longitudinal data to data with one variable per subject
vars_static <- setdiff(vars_of_interest, c("crp", "pct"))
data.tind <- data.full %>% 
  select(all_of(c(metrange, vars_of_interest, "subject.id", "day"))) %>% 
  filter(day %in% c(0, 30)) %>% 
  pivot_longer(-c(all_of(vars_static), subject.id, day), names_to = "metabolite") %>%
  pivot_wider(id_cols = c(subject.id, metabolite, all_of(vars_static)), names_from = day) %>% 
  mutate(met_d30_d0 = `30` - `0`) %>% 
  pivot_wider(names_from = metabolite, values_from = met_d30_d0, id_cols = c(subject.id, all_of(vars_static)))

#Fit a linear model with for each metabolite
results_lm <- NULL
for (i in 1:length(metrange)) {
  dat_lm <- data.tind %>% 
    select(all_of(c(metrange[i], vars_of_interest)))
  
  l_model <- lm(data = dat_lm, as.formula(paste(metrange[i], "~", paste(vars_of_interest, collapse = " + "))))
  results_lm <- rbind.data.frame(results_lm, data.frame(metabolite = metrange[i], summary(l_model)$coefficients[-1, ]) %>% rownames_to_column("covariate"))
}




#get top 10s out
results_lm %>%
  filter(covariate == "crp") %>% 
  dplyr::arrange(desc(abs(t.value))) %>% slice_head(n = 10)



#hospitalization time analysis

#linear model with day1-0, day2-0, day4-0, day30-0
results_lm_hosp <- NULL
vars_hosp <- c("age", "sex", paste0("day", c(1, 2, 4, 30), "_0"))
for (i in 1:length(metrange)) {
  dat_lm <- data.full %>% 
    select(all_of(c(metrange[i], "age", "sex", "day", "subject.id", "hospitalization.time"))) %>% 
    #feature engineering
    pivot_wider(names_from = day, values_from = all_of(metrange[i]), names_prefix = "day") %>% 
    mutate(day1_0 = day1 - day0, day2_0 = day2 - day0, day4_0 = day4 - day0, day30_0 = day30 - day0)
  
  lm_model_hosp <- lm(data = dat_lm, as.formula(paste("hospitalization.time ~", paste(c("age", "sex", paste0("day", c(1, 2, 4, 30), "_0")), collapse = " + "))))
  results_lm_hosp <- rbind.data.frame(results_lm_hosp, data.frame(metabolite = metrange[i], summary(lm_model_hosp)$coefficients[-1, ]) %>% rownames_to_column("covariate"))
}


save(results_lm, results_lmm, results_lm_hosp, file = "results/linear_model_results.Rdata")


#load results from script:
load("results/linear_model_results.Rdata")

#ordinal regression model with day1-0, day2-0, day4-0, day30-0 on curb
results_lm_curb <- NULL
vars_curb <- c("age", "sex", paste0("day", c(1, 2, 4, 30), "_0"))
for (i in 1:length(metrange)) {
  dat_lm <- data.full %>% 
    select(all_of(c(metrange[i], "age", "sex", "day", "subject.id", "curb"))) %>% 
    #feature engineering
    pivot_wider(names_from = day, values_from = all_of(metrange[i]), names_prefix = "day") %>% 
    mutate(day1_0 = day1 - day0, day2_0 = day2 - day0, day4_0 = day4 - day0, day30_0 = day30 - day0) %>% 
    mutate(curb = as.factor(curb), sex = as.factor(sex))
  
  lm_model_curb <- MASS::polr(data = dat_lm, as.formula(paste("curb ~", paste(c("age", paste0("day", c(1, 2, 4, 30), "_0")), collapse = " + "))))
  results_lm_curb <- rbind.data.frame(results_lm_curb, data.frame(metabolite = metrange[i], summary(lm_model_hosp)$coefficients[-1, ]) %>% rownames_to_column("covariate"))
}


lm_model_curb$
MASS::polr()
test_om <- MASS::polr(data = dat_lm,  as.formula(paste("curb ~", paste(c("age", paste0("day", c(1, 2, 4, 30), "_0")), collapse = " + "))), Hess = TRUE)
coef(summary(test_om))
test_om$edf
confint(profile(test_om))
pairs(profile(test_om))

11.6227210 + 2*qt(0.025, 1.1342271, 12)
11.6227210 + 7.882866
