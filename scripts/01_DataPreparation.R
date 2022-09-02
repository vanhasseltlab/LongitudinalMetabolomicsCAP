# Data preparation

# Load libraries
library(tidyverse)

# Load functions
source("functions/AllFunctions.R")

# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")

# Data preparation
#Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)

#Remove metabolites with to many NAs
data.clean <- DataCleaning(data.reduced, attr(data.reduced, "metrange"))

#Add longitudinal markers: CRP and PCT
data.full <- AddCRPAndPCT(data.reduced)

#Add ratios and sums to data.full
data.full.clean <- AddRatiosAndSums(AddCRPAndPCT(data.clean))

#Log transform and scale metabolite data
data.pretreated <- DataPretreatment(data.full.clean, scaling = "auto", 
                                    metrange = c(attr(data.full.clean, "metrange"), attr(data.full.clean, "ratiorange"))) %>% 
  mutate(subject.id_orig = subject.id, subject.id = paste("Patient", as.numeric(as.factor(subject.id)))) %>% 
  mutate(subject.id = factor(subject.id, levels = paste("Patient", 1:length(unique(data.full.clean$subject.id)))))

#Save data
#save(data.full.clean, data.pretreated, file = "data/01_data_clean.Rdata")
save(data.pretreated, file = "data/01_data_clean.Rdata")
