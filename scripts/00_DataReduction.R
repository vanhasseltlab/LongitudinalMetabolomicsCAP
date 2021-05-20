## Code for data reduction: 
##  Removing unnecessary variables and translating text to english

## Setup R --------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Load libraries
library(plyr)
library(dplyr)

## Input data -----------------------------------------------------------------
# Load metabolomics data after quality control
data.raw <- read.csv2("data/00_data_raw.csv")

## Data reduction--------------------------------------------------------------
# Reduce dataset, only keep nesseccary variables
metabolomicsdata <- data.raw[, 6:374]
patientdata <- select(data.raw, 
                      "sample.id"= NMC.name,
                      "subject.id" = Studynr,
                      "day" = Day,
                      "pathogen" = verwekker_1,
                      "pathogen.group" = verwekker_groep, 
                      "age" = leeftijd,  
                      "sex" = geslacht,
                      "psi.score" = Fine_klasse,
                      "nursing.home.resident" = VPH, # verpleeghuisopname
                      "renal.disease" = nierziekte,
                      "liver.disease" = leverziekte,
                      "congestive.heart.failure" = hartfalen,
                      "cns.disease" = CNS,
                      "malignancy" = maligniteit,
                      "altered.mental.status" = verward, 
                      "respiratory.rate" = Afreq,
                      "systolic.blood.pressure" = RR_syst,
                      "temperature" = temp,
                      "pulse" = pols,
                      "pH" = pH_0,
                      "BUN" = ureum_0, # BUN = Blood Urea Nitrogen concentration
                      "sodium" = natrium_0, 
                      "glucose" = glucose_0,
                      "hematocrit" = Ht_0,
                      "partial.pressure.of.oxygen" = PO2_0,
                      "pleural.effusion.on.x-ray" = Xth_pv,
                      "race" = Ras, 
                      "duration.of.symptoms.before.admission" = duur_sympt,
                      "antibiotic.treatment.before.admission" = AB_thuis,
                      "corticosteroid.use.before.admission" = med_cort,
                      "COPD" = COPD, 
                      "diabetes" = DM,
                      "oxygen.saturation" = sat,
                      "supplemental.oxygen.required" = O2,
                      "antibiotics.started.in.hospital" = gestart_AB_ZH, 
                      "CRP" = CRP_0, 
                      "leukocyte count" = leuko_0,
                      "IL-6" = IL6_0,
                      "IL-10" = IL10_0)

data.reduced <- bind_cols(metabolomicsdata, patientdata)


data.reduced$pathogen.group <- revalue(data.reduced$pathogen.group, c(
  "Streptococcus pneumoniae" = "S. pneumoniae"))

data.reduced$psi.score <- revalue(data.reduced$psi.score, c(
  "0" = "0-50",
  "<70" = "51-70", 
  "71-90" = "71-90", 
  "91-130" = "91-130", 
  ">130" = "131-395"))
data.reduced$psi.score <- factor(data.reduced$psi.score, levels = c("0-50", "51-70", "71-90","91-130", "131-395"))

#Translate Dutch to English
data.reduced$sex <- revalue(data.reduced$sex, c("man" = "Male", "vrouw" = "Female"))
data.reduced$race <- revalue(data.reduced$race, c("wit" = "White"))

data.reduced <- data.frame(lapply(data.reduced, function(x) { gsub("nee", "No", x) }))
data.reduced <- data.frame(lapply(data.reduced, function(x) { gsub("ja", "Yes", x) }))


# Define metabolite and metadata range in data after quality control
metadata <- c("sample.id", "subject.id", "day", "pathogen", "pathogen.group", "age", "sex", "psi.score", "nursing.home.resident",
              "renal.disease", "liver.disease", "congestive.heart.failure", "cns.disease", "malignancy", "altered.mental.status",
              "respiratory.rate", "systolic.blood.pressure", "temperature", "pulse", "pH" , "BUN" , "sodium", "glucose", 
              "hematocrit", "partial.pressure.of.oxygen", "pleural.effusion.on.x.ray", "race", "duration.of.symptoms.before.admission",
              "antibiotic.treatment.before.admission", "corticosteroid.use.before.admission", "COPD", "diabetes", 
              "oxygen.saturation", "supplemental.oxygen.required", "antibiotics.started.in.hospital", "CRP", "leukocyte.count",
              "IL.6", "IL.10")

psi.score.components <- c("nursing.home.resident",                "renal.disease",                         "liver.disease",                        
                          "congestive.heart.failure",              "cns.disease",                           "malignancy",                           
                          "altered.mental.status",                 "respiratory.rate",                      "systolic.blood.pressure",              
                          "temperature",                           "pulse",                                 "pH",                                   
                          "BUN",                                   "sodium",                               "glucose",                              
                          "hematocrit",                            "partial.pressure.of.oxygen",            "pleural.effusion.on.x.ray"  )

infection_markers <- c("CRP", "leukocyte.count", "IL.6", "IL.10")

covnames <- c("age", "sex", psi.score.components, "race", "duration.of.symptoms.before.admission", 
              "antibiotic.treatment.before.admission", "corticosteroid.use.before.admission", "COPD", "diabetes")

metrange <- setdiff(names(data.reduced), metadata)

# Define metabolite data as numeric data for analysis
data.reduced[,metrange] <- apply(data.reduced[,metrange], 2, function(x) {
  as.numeric(as.character(x))
})

## Output data ----------------------------------------------------------------

# Save reduced metabolomics dataframe in .csv
write.csv2(data.reduced, file = "data/00_data_reduced.csv", row.names = FALSE)



