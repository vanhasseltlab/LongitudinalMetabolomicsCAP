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

#Select only relevant columns
all_used_columns <- c("day", "subject.id", "hospitalization.time", "curb", "pct", "crp", "crea",
                      attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))
data.pretreated <- data.pretreated %>% select(all_of(all_used_columns))

#Save data
#save(data.full.clean, data.pretreated, file = "data/01_data_clean.Rdata")
save(data.pretreated, file = "data/01_data_clean.Rdata")
write.table(data.pretreated, file = "data/01_data_clean.csv", row.names = FALSE, quote = T, sep = ",")


#Create table for patient characteristics
## Patient characteristics ---
patient_chars <- c("age", "sex", "curb", "renal.disease", "congestive.heart.failure",
                   "malignancy", "COPD", "diabetes", "duration.of.symptoms.before.admission",
                   "antibiotic.treatment.before.admission", "corticosteroid.use.before.admission", 
                   "hospitalization.time")
patient_chars_labels <-  c("Age (years)", "Sex", "CURB score", "Kidney disease", "Cardiovascular disease",
                           "Malignancy", "COPD", "Diabetes", "Duration of symptoms before admission (days)", 
                           "Antibiotic treatment before admission", "Corticosteroid use before admission", 
                           "Length of stay (days)")
data_one_obs <- data.full.clean %>% 
  distinct(subject.id, .keep_all = TRUE) %>% 
  select(all_of(patient_chars), pathogen.group)

variables <- as.list(patient_chars_labels)
names(variables) <- patient_chars
labels_list <- list(variables = variables,
                    groups = list("S. pneumoniae"))
table1::table1(list(`CAP patients` = data_one_obs), labels = labels_list)
