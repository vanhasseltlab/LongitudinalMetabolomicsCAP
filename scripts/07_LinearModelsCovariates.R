#Regression of metabolites on different covariates

# Load libraries
library(tidyverse)
library(lme4)
# Source functions
source("functions/AllFunctions.R")
source("functions/DataCleaning.R")

#Data
# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
#Add longitudinal markers: CRP and PCT
data.full <- AddCRPAndPCT(data.reduced)

#Covariates of interest:
#age, sex, antibiotica bij start, antibiotica switch, CURB, CRP, PCT, 
# duration of symptoms before admission, hospitalization time
#Note: "duration.of.symptoms.before.admission", "antibiotic.treatment.before.admission" have a lot of missingness

#Select variables of interest
vars_of_interest <- c("age", "sex", "curb", "hospitalization.time", "crp", "pct")

#Define metrange
metrange <- attr(data.full, "metrange")

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
