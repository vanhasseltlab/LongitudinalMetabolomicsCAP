## This is code to run all analysis for the Longitunal Metabolomics Project
## Code developped by Ilona den Hartog and Laura Zwep

# Code contents
# - Data cleaning
# - Data pretreatment
# - Data analysis for the comparisons:
#   - Comparison 1: S. pneumoniae versus atypical bacteria + virusses
#   - Comparison 2: Atypical bacteria versus S. pneumoniae + virusses
#   - Comparison 3: Virusses versus S. pneumoniae + atypical bacteria
# - Statistical data analysis:
#   - T-test + FDR multiple testing correction for comparison of 1 pathogen 
#      group versus the remaining 2 pathogen groups. 
  

## Setup R --------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Load libraries
# library(tidyverse)
library(dplyr)
# library(ggplot2)

# Source functions
source("functions/DataCleaning.R")
# source("functions/LeidenColoring.R")


## Input data -----------------------------------------------------------------
dat.raw <- read.csv2("data/00_data_raw.csv")
# Load metabolomics data after quality control
data.reduced <- read.csv2("data/00_data_reduced.csv")

# Define metabolite and metadata range in data after quality control
nonmetrange <- c("sample.id", "subject.id", "day", "pathogen", "pathogen.group", "age", "sex", "psi.score", "nursing.home.resident",
                 "renal.disease", "liver.disease", "congestive.heart.failure", "cns.disease", "malignancy", "altered.mental.status",
                 "respiratory.rate", "systolic.blood.pressure", "temperature", "pulse", "pH" , "BUN" , "sodium", "glucose", 
                 "hematocrit", "partial.pressure.of.oxygen", "pleural.effusion.on.x.ray", "race", "duration.of.symptoms.before.admission",
                 "antibiotic.treatment.before.admission", "corticosteroid.use.before.admission", "COPD", "diabetes", 
                 "oxygen.saturation", "supplemental.oxygen.required", "antibiotics.started.in.hospital", "CRP", "leukocyte.count",
                 "IL.6", "IL.10")
metrange <- setdiff(names(data.reduced), nonmetrange)

infection_markers <- c("CRP", "leukocyte.count", "IL.6", "IL.10")

psi.score.components <- c("nursing.home.resident",                "renal.disease",                                                 
                          "congestive.heart.failure",              "cns.disease",                           "malignancy",                           
                          "altered.mental.status",                 "respiratory.rate",                      "systolic.blood.pressure",              
                          "temperature",                           "pulse",                                 "pH",                                   
                          "BUN",                                   "sodium",                               "glucose",                              
                          "hematocrit",                            "partial.pressure.of.oxygen",            "pleural.effusion.on.x.ray"  )

# cov <- c("age", "sex", psi.score.components, "duration.of.symptoms.before.admission", 
#               "antibiotic.treatment.before.admission", "COPD", "diabetes")

# Save colnames
save(nonmetrange, metrange, infection_markers, psi.score.components, file = "data/column_names.Rdata")

## Data cleaning --------------------------------------------------------------
data.clean <- DataCleaning(data.reduced, metrange, nonmetrange)
#Reset metrange after cleaning
metrange <- setdiff(names(data.clean), nonmetrange)

## Data pretreatment ----------------------------------------------------------
# Log2 transformation of metabolite values +1. 
logdata <- data.clean
logdata[, metrange] <- apply(logdata[, metrange], MARGIN = 2, FUN = function(x) { 
  log2(x + 1)
})

# Autoscaling of log transformed data. 
data.pretreated <- logdata
data.pretreated[, metrange] <- apply(data.pretreated[, metrange], MARGIN = 2, FUN = function(x) {
  (x - mean(x)) / sd(x)
})

## Data summary --------------------------------------------------------------
dat <- data.pretreated

View(select(dat, subject, day))


## Data imputation -----------------------------------------------------------



