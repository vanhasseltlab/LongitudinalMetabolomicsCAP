# Prepare data by combining patient data with measurements form different metabolomic platforms 

# Load libraries
library(foreign)
library(readr)
library(dplyr)
library(tidyr)

## Input data -----------------------------------------------------------------

# Patient data from Ovidius and 3P studies
patdataraw <- read.spss("data/raw/np_171031_SV_db02_complete_studydatabase_Ovidius_3P_CAP_studies.sav", 
                        use.value.labels = TRUE, to.data.frame = TRUE)

patdataraw_extra <- read.spss("data/raw/Aanvullende data illona SPSS.sav", 
                              use.value.labels = TRUE, to.data.frame = TRUE)
# Remove duplicate row (Ovidius 115, row 28)
patdataraw_extra <- patdataraw_extra[-28,]

# Merge patient data
patdataraw_merged <- left_join(patdataraw, patdataraw_extra, by = c("studie" = "Studie", "Studienr"))

# Metabolomics data
# - data with and without caution
# - Manual correction of one incomplete sample number (s-...)

# Amines
nmcdatarawAmines <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_amines.csv", 
                               ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                               trim_ws = TRUE, skip = 3)
nmcdatarawAmines[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAminesCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_amines_caution.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 3)
nmcdatarawAminesCaut[190,1] <- "303, S341185, 17-09-2010" 

# Acyl carnitines
nmcdatarawAcyl <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_acylcarnitines.csv", 
                             ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                             trim_ws = TRUE, skip = 3) %>%
  select(-2)
nmcdatarawAcyl[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAcylCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_acylcarnitines_caution.csv", 
                                 ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                 trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcylCaut[190,1] <- "303, S341185, 17-09-2010" 

# Organic Acids
nmcdatarawAcids <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_organic_acids.csv", 
                              ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                              trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcids[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawAcidsCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_organic_acids_caution.csv", 
                                  ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawAcidsCaut[190,1] <- "303, S341185, 17-09-2010" 

# Negative Lipids
nmcdatarawNegLip <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_negative_lipids.csv", 
                               ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                               trim_ws = TRUE, skip = 5)%>% 
  select(-2)
nmcdatarawNegLip[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawNegLip)[names(nmcdatarawNegLip) == "Original name"] <- "Original"

nmcdatarawNegLipCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_negative_lipids_caution.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 5)%>% 
  select(-2)
nmcdatarawNegLipCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawNegLipCaut)[names(nmcdatarawNegLipCaut) == "Original name"] <- "Original"

# Positive Lipids
# TG's
nmcdatarawPosLipTG <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_tg.csv", 
                                 ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                 trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipTG[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipTG)[names(nmcdatarawPosLipTG) == "sampleID"] <- "Original"

nmcdatarawPosLipTGCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_tg_caution.csv", 
                                     ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                     trim_ws = TRUE, skip = 4)%>% 
  select(-2)
nmcdatarawPosLipTGCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipTGCaut)[names(nmcdatarawPosLipTGCaut) == "sampleID"] <- "Original"
# non TG's 
nmcdatarawPosLipNonTG <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_non-tg.csv", 
                                    ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                    trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipNonTG[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipNonTG)[names(nmcdatarawPosLipNonTG) == "sampleID"] <- "Original"

nmcdatarawPosLipNonTGCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_positive_lipids_non-tg_caution.csv", 
                                        ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                        trim_ws = TRUE, skip = 4)%>% 
  select(-c(2,3))
nmcdatarawPosLipNonTGCaut[187,1] <- "303, S341185, 17-09-2010" 
names(nmcdatarawPosLipNonTGCaut)[names(nmcdatarawPosLipNonTGCaut) == "sampleID"] <- "Original"

# Signalling molecules (oxidative stress & oxylipins)
# High pH
nmcdatarawSignHighPH <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_high_pH.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE, skip = 3)%>% 
  select(-2, -length(colnames(.))) # remove also last column because it is empty
nmcdatarawSignHighPH[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawSignHighPHCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_high_pH_caution.csv", 
                                       ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                       trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignHighPHCaut[190,1] <- "303, S341185, 17-09-2010" 

# Low PH
nmcdatarawSignLowPH <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_low_pH.csv", 
                                  ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignLowPH[190,1] <- "303, S341185, 17-09-2010" 

nmcdatarawSignLowPHCaut <- read_delim("data/raw/NMC-17-058_ZonMW-Ilona_Data_v3(1)_signalling_low_pH_caution.csv",
                                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                      trim_ws = TRUE, skip = 3)%>% 
  select(-2)
nmcdatarawSignLowPHCaut[190,1] <- "303, S341185, 17-09-2010" 

# 1.c) Upload connecting files: 
# Semi-randomized sample list containing sampleID, study and study number and day after admission
connectraw <- read_csv2("data/raw/20180405_IDH_semi-randomized_sample_list_V4.csv")
linkfile <- subset(connectraw, select=c(Study_cohort, Studynr, Day, SampleID))
linkfile$SampleID <- as.character(linkfile$SampleID)
linkfile[193,4] <- "303, S341185, 17-09-2010"

# Batch design
batchDesign <- read_delim("data/raw/NMC-17-058_Batch_Design.csv", 
                          ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                          trim_ws = TRUE) %>% 
  select(-3, -4, 6)
batchDesign[192,2] <- "303, S341185, 17-09-2010" 


## Data assembly --------------------------------------------------------------
# All data
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
  left_join(patdataraw_merged,           by=c("Study_cohort" = "studie", "Studynr" = "Studienr"))  

# Remove duplicate metabolites (duplicate metabolite in NegLips and double measured metabolites because of extended signalling lipids platform)
alldata_corrected <- alldata %>% 
  select(-c(`NMC name.y`, `FA (17:0).y`, AEA.y, LEA.y, O.AEA.y, SEA.y, DEA.y)) %>% 
  rename(`NMC name` = `NMC name.x`,
         `FA (17:0)` = `FA (17:0).x`, 
         AEA = AEA.x, 
         LEA = LEA.x, 
         O.AEA = O.AEA.x, 
         SEA = SEA.x, 
         DEA = DEA.x)

# Remove sample causative pathogen H. Influenzae (H. Influenzae wrongly noted as viral pathogen in patient data)
alldata_corrected <- filter(alldata_corrected, verwekker_1 != "H. Influenza")

# # Diagnosisdata 
# diagnosisdata <- subset(alldata_corrected, Day==0) 

# Treatment resonse data (longitudinal) and rename IDs
longitudinaldata <- alldata_corrected %>%
  group_by(Studynr, Study_cohort) %>%
  filter(length(Study_cohort) > 1) %>%
  ungroup()

## Output data ----------------------------------------------------------------
# All data
write.csv2(alldata_corrected, file = "data/00_data_raw_all.csv", row.names = FALSE)
#Longitudinal data
write.csv2(longitudinaldata, file = "data/00_data_raw.csv", row.names = FALSE)


  