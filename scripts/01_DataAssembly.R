## Upload raw metabolomics data, couple to patient data, and clean dataset

## Setup R -----
# Version information:
# * R version 3.4.2
# * R-studio version 1.0.153

# Load libraries
library(foreign)
library(readr)
library(dplyr)
library(tidyr)

# Change working directory
setwd("~/LongitudinalMetabolomics/01_DataAssembly")

# Clean global environment
rm(list=ls())

## Input data----

## Data folder pathway 
datafolder = "~/ZonMW/5. Raw data/"

# Upload patient data (.sav file)
patdatasav <- paste(datafolder, "np_171031_SV_db02_complete_studydatabase_Ovidius_3P_CAP_studies.sav", sep="")
patdataraw <- read.spss(patdatasav, use.value.labels=TRUE, to.data.frame=TRUE)

# Upload metabolomics data (nmc data)
# - Metabolomics data and metabolomics data with caution
# - Manually corrected incomplete s-number of 1 sample in order to allow leftjoin by SampleID
# - Manually removed columns that were interfering with leftjoin (different NMC names for different platforms) or not used (in BatchDesign)

# Load metabolomics data after quality control
# Amines
nmcdatarawAmines <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_amines.csv", sep=""), 
                               ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                               trim_ws = TRUE, skip = 3)
nmcdatarawAmines[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAminesCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_amines_caution.csv", sep=""), 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 3)
nmcdatarawAminesCaut[190,1] <- "303, S341185, 17-09-2010" 

# Acyl carnitines
nmcdatarawAcyl <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_acylcarnitines.csv", sep=""), 
                             ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                             trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcyl[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAcylCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_acylcarnitines_caution.csv", sep=""), 
                                 ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                 trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcylCaut[190,1] <- "303, S341185, 17-09-2010" 

# Organic Acids
nmcdatarawAcids <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_organic_acids.csv", sep=""), 
                              ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                              trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcids[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAcidsCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_organic_acids_caution.csv", sep=""), 
                                  ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcidsCaut[190,1] <- "303, S341185, 17-09-2010" 

# Negative Lipids
nmcdatarawNegLip <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_negative_lipids.csv", sep=""), 
                               ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                               trim_ws = TRUE, skip = 5)%>% 
  select(-2)
nmcdatarawNegLip[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawNegLip)[names(nmcdatarawNegLip) == "Original name"] <- "Original"

nmcdatarawNegLipCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_negative_lipids_caution.csv", sep=""), 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 5)%>% 
  select(-2)
nmcdatarawNegLipCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawNegLipCaut)[names(nmcdatarawNegLipCaut) == "Original name"] <- "Original"

# Positive Lipids
# TG's
nmcdatarawPosLipTG <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_tg.csv", sep=""), 
                                 ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                 trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipTG[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipTG)[names(nmcdatarawPosLipTG) == "sampleID"] <- "Original"

nmcdatarawPosLipTGCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_tg_caution.csv", sep=""), 
                                     ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                     trim_ws = TRUE, skip = 4)%>% 
  select(-2)
nmcdatarawPosLipTGCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipTGCaut)[names(nmcdatarawPosLipTGCaut) == "sampleID"] <- "Original"
# non TG's 
nmcdatarawPosLipNonTG <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_non-tg.csv", sep=""), 
                                    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                    trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipNonTG[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipNonTG)[names(nmcdatarawPosLipNonTG) == "sampleID"] <- "Original"

nmcdatarawPosLipNonTGCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_non-tg_caution.csv", sep=""), 
                                        ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                        trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipNonTGCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipNonTGCaut)[names(nmcdatarawPosLipNonTGCaut) == "sampleID"] <- "Original"

# Signalling molecules (oxidative stress & oxylipins)
# High pH
nmcdatarawSignHighPH <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_high_pH.csv", sep=""), 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 3)%>% 
  select(-2, -length(colnames(.))) # remove also last column because it is empty
nmcdatarawSignHighPH[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawSignHighPHCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_high_pH_caution.csv", sep=""), 
                                       ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                       trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignHighPHCaut[190,1] <- "303, S341185, 17-09-2010" 

# Low PH
nmcdatarawSignLowPH <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_low_pH.csv", sep=""), 
                                  ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignLowPH[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawSignLowPHCaut <- read_delim(paste(datafolder, "NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_low_pH_caution.csv", sep=""), 
                                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                      trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignLowPHCaut[190,1] <- "303, S341185, 17-09-2010" 

# Upload connecting files
# Semi-randomized sample list containing sampleID, study and study number and day after admission
connectraw <- read.csv2(paste(datafolder,"20180405_IDH_semi-randomized sample list_V4.csv", sep=""), fileEncoding = "latin1")
linkfile <- subset(connectraw, select=c(Study.1, Studynr, Day, SampleID))
linkfile$SampleID <- as.character(linkfile$SampleID)
linkfile[193,4] <- "303, S341185, 17-09-2010"

# Batch design containing batch information
batchDesign <- read_delim(paste(datafolder,"NMC-17-058_Batch_Design.csv", sep=""), 
                          ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                          trim_ws = TRUE) %>% 
  select(-3, -4, 6)
batchDesign[192,2] <- "303, S341185, 17-09-2010" 

## Generate dataframe containing metabolomics and patient data ----
# All data (alldata)
joinby <- c("SampleID" = "Original")
alldata <- left_join(linkfile, nmcdatarawAmines, by=joinby) %>% 
  left_join(nmcdatarawAminesCaut,     by=joinby) %>% 
  left_join(nmcdatarawAcyl,           by=joinby) %>% 
  left_join(nmcdatarawAcylCaut,       by=joinby) %>% 
  left_join(nmcdatarawAcids,          by=joinby) %>% 
  left_join(nmcdatarawAcidsCaut,      by=joinby) %>% 
  left_join(nmcdatarawNegLip,         by=joinby) %>% 
  left_join(nmcdatarawNegLipCaut,     by=joinby) %>% 
  left_join(nmcdatarawSignHighPH,     by=joinby) %>% 
  left_join(nmcdatarawSignHighPHCaut, by=joinby) %>%   
  left_join(nmcdatarawSignLowPH,      by=joinby) %>% 
  left_join(nmcdatarawSignLowPHCaut,  by=joinby) %>% 
  left_join(nmcdatarawPosLipTG,       by=joinby) %>% 
  left_join(nmcdatarawPosLipTGCaut,   by=joinby) %>% 
  left_join(nmcdatarawPosLipNonTG,    by=joinby) %>% 
  left_join(nmcdatarawPosLipNonTGCaut,by=joinby) %>% 
  left_join(batchDesign,          by=c("SampleID"="sampleID")) %>% 
  left_join(patdataraw,           by=c("Study.1" = "studie", "Studynr" = "Studienr"))  

# Remove duplicate metabolites (duplicate metabolite in NegLips and double measured metabolites because of extended signalling lipids platform)
alldata_corrected_1 <- alldata %>% 
  select(-c(`NMC name.y`, `FA (17:0).y`, AEA.y, LEA.y, O.AEA.y, SEA.y, DEA.y)) %>% 
  rename(`NMC name` = `NMC name.x`,
         `FA (17:0)` = `FA (17:0).x`, 
         AEA = AEA.x, 
         LEA = LEA.x, 
         O.AEA = O.AEA.x, 
         SEA = SEA.x, 
         DEA = DEA.x)

# Remove sample causative pathogen H. Influenzae (H. Influenzae wrongly noted as viral pathogen in patient data)
alldata_corrected <- filter(alldata_corrected_1, verwekker_1 != "H. Influenza")

# Information about which columns contain metabolites and of which platforms (Note: this is manually looked up and added)
# Total of 374-5 = 369 metabolites measured :)
metrange1    <- 6:374 #start column metabolites : final column containing metabolites 
# Structure of the data
samples      <- 1:5      # 5 columns
#metrange
patdat       <- 381:436  # 55 columns

# Diagnosisdata 
diagnosisdata_afterQC <- subset(alldata_corrected, Day==0) 

# Treatment resonse data (timecourse)
timecoursedata_afterQC <- alldata_corrected %>% 
  group_by(Studynr,Study.1) %>% 
  filter(length(Study.1)>1) %>%
  ungroup()  

## Output data ----
#Alldata
write.csv2(alldata_corrected, file = paste(Sys.Date(), "_IDH_alldata_afterQC.csv", sep=""), row.names = FALSE)
save(alldata_corrected, file = paste(Sys.Date(), "_IDH_alldata_afterQC.Rdata", sep=""))
#Diagnosisdata
write.csv2(diagnosisdata_afterQC, file = paste(Sys.Date(), "_IDH_diagnosisdata_afterQC.csv", sep=""), row.names = FALSE)
save(diagnosisdata_afterQC, file = paste(Sys.Date(), "_IDH_diagnosisdata_afterQC.Rdata", sep=""))
#Timecoursedata
write.csv2(timecoursedata_afterQC, file = paste(Sys.Date(), "_IDH_timecoursedata_afterQC.csv", sep=""), row.names = FALSE)
save(timecoursedata_afterQC, file = paste(Sys.Date(), "_IDH_timecoursedata_afterQC.Rdata", sep=""))
