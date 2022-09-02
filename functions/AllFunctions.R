## This script contains all functions nessecary to reproduce the data from the Longitudinal Metabolomics project
library(foreign)
library(openxlsx)


#### Remove redundant variables from raw data using reduction function ----
ReduceData <- function(data.raw){
  metabolomicsdata <- data.raw[, 6:374]
  patient_vars <- c("sample.id" = "NMC.name", # Metabolomics sample id
                    "sample.id.S" = "SampleID", # Hospital sample id
                    "subject.id" = "Studynr", # Study patient id
                    "day" = "Day",
                    "pathogen" = "verwekker_1",
                    "pathogen.group" = "verwekker_groep", 
                    "age" = "leeftijd",  
                    "sex" = "geslacht",
                    "psi.score" = "Fine_klasse",
                    "nursing.home.resident" = "VPH", # verpleeghuisopname
                    "renal.disease" = "nierziekte",
                    "liver.disease" = "leverziekte",
                    "congestive.heart.failure" = "hartfalen",
                    "cns.disease" = "CNS",
                    "malignancy" = "maligniteit",
                    "altered.mental.status" = "verward", 
                    "respiratory.rate" = "Afreq",
                    "systolic.blood.pressure" = "RR_syst",
                    "temperature" = "temp",
                    "pulse" = "pols",
                    "pH" = "pH_0",
                    "BUN" = "ureum_0", # BUN = Blood Urea Nitrogen concentration
                    "sodium" = "natrium_0", 
                    "glucose" = "glucose_0",
                    "hematocrit" = "Ht_0",
                    "partial.pressure.of.oxygen" = "PO2_0",
                    "pleural.effusion.on.x-ray" = "Xth_pv",
                    "race" = "Ras", 
                    "duration.of.symptoms.before.admission" = "duur_sympt",
                    "antibiotic.treatment.before.admission" = "AB_thuis",
                    "corticosteroid.use.before.admission" = "med_cort",
                    "COPD" = "COPD", 
                    "diabetes" = "DM",
                    "oxygen.saturation" = "sat",
                    "supplemental.oxygen.required" = "O2",
                    "antibiotics.started.in.hospital" = "gestart_AB_ZH", 
                    "CRP" = "CRP_0", 
                    "leukocyte count" = "leuko_0",
                    "IL-6" = "IL6_0",
                    "IL-10" = "IL10_0",
                    "hospitalization.time" = "Opnameduur",
                    "death" = "dood_ontslag")
  
  patientdata <- select(data.raw, all_of(patient_vars))
  
  data.reduced <- bind_cols(metabolomicsdata, patientdata)
  
  data.reduced$pathogen.group <- recode(data.reduced$pathogen.group, `Streptococcus pneumoniae` = "S. pneumoniae")
  
  data.reduced$psi.score <- recode(data.reduced$psi.score, 
                                   `0` = "0-50", `<70` = "51-70", 
    `71-90` = "71-90", 
    `91-130 `= "91-130", 
    `>130` = "131-395")
  data.reduced$psi.score <- factor(data.reduced$psi.score, levels = c("0-50", "51-70", "71-90","91-130", "131-395"))
  
  #Translate Dutch to English
  data.reduced$sex <- recode(data.reduced$sex, man = "Male", vrouw = "Female")
  data.reduced$race <- recode(data.reduced$race, wit = "White")
  data.reduced$death <- recode(data.reduced$death, levend = "alive", dood = "dead")
  
  ja_nee_columns <- which(apply(data.reduced, 2, function(x) any(x == "ja" | x == "nee", na.rm = T)))
  
  data.reduced[, ja_nee_columns] <- data.frame(lapply(data.reduced[, ja_nee_columns], function(x) { 
    x <- gsub("nee", "No", x)
    x <- gsub("ja", "Yes", x) 
    return(x)
  }))
  
  ## Calculate CURB score at day 0 for all patients (subject.id's) and add to reduced data
  data.reduced$curb <- with(data.reduced, {
    (BUN > 7) + (as.character(altered.mental.status) == "Yes") +
    (respiratory.rate >= 30) + (systolic.blood.pressure < 90)})
  
  #make subject.id character
  data.reduced$subject.id <- as.character(data.reduced$subject.id)
  
  # Define metabolite and metadata range in data after quality control
  metrange <- setdiff(names(data.reduced), c(names(patient_vars), "curb"))
  attr(data.reduced, "metrange") <- metrange
  
  return(data.reduced)
}

DataCleaning <- function(dat, metrange) {
  # This function removes metabolites (columns) with missing values. And patients with
  # death outcome
  nonmetrange <-  setdiff(names(dat), metrange)
  
  dat1 <- dat[dat$death == "alive", ]
  # Subset metabolites from partly cleaned dataset
  dat.subset.2 <- dat1[, metrange]
  # Calculate the numer of NA's per column
  nacols <- apply(dat.subset.2, 2, function(X) {sum(is.na(X))})
  # keep column names of metabolites without missing metabolites
  keep <- colnames(dat.subset.2[, nacols < 5])
  
  # Combine cleaned metabolomicsdat with non-metabolomicsdat 
  #  and return clean dat set.
  dat.clean <- cbind(dat1[, keep], dat1[, nonmetrange])
  attr(dat.clean, "metrange") <- keep
  
  return(dat.clean)
}

DataPretreatment <- function(dat, metrange, scaling = "auto") {
  # This function applies log transformation and auto or pareto scaling to metabolite data
  nonmetrange <-  setdiff(names(dat), metrange)
  # Log2 transformation of metabolite values +1. 
  logdata <- dat
  logdata[, metrange] <- apply(logdata[, metrange], MARGIN = 2, FUN = function(x) { 
    log2(x + 1)
  })
  if (scaling == "auto"){ 
    # Autoscaling of log transformed data. 
    data.pretreated <- logdata
    data.pretreated[, metrange] <- apply(data.pretreated[, metrange], MARGIN = 2, FUN = function(x) {
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    })
  } else if (scaling == "pareto") {
    # Paretoscaling of log transformed data (using the square root of the SD)
    data.pretreated <- logdata
    data.pretreated[, metrange] <- apply(data.pretreated[, metrange], MARGIN = 2, FUN = function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    })
  } else { 
    print("Please set scaling to auto or pareto")}
  return(data.pretreated)
}

#### Function to calculate PCT levels from CLIA assays ----
## Function to calculate PCT values from plate reader (pr) results
PCTcalculator <- function(pr_results, pm_location){
  #Create return list and 
  return_list <- list()
  
  # Adjust formatting of plate reader results for further analysis
  colnames(pr_results) <- c("location", "RLU")
  pr_results <- separate(pr_results, col = location, into = c("RowWell", "ColWell"), sep = 1)
  pr_results$ColWell <- as.character(as.numeric(pr_results$ColWell))
  
  pr_results_mapped <- pr_results %>% 
    pivot_wider(names_from = ColWell, values_from = RLU) %>% 
    column_to_rownames(var = "RowWell") 
  
  # Adjust formatting of plate map to plate reader results
  pm_sample_ids <- read.xlsx(pm_location, rows = c(9:17), cols = c(1:13), 
                             colNames = T, rowNames = T) %>%   
    rownames_to_column("RowWell") %>% 
    pivot_longer(-RowWell, names_to = "ColWell", values_to = "sample_id", values_transform = list(sample_id=as.character)) 
  
  pm_concentrations <- read.xlsx(pm_location, rows = c(20:28), cols = c(1:13), 
                                 colNames = T, rowNames = T) %>% 
    rownames_to_column("RowWell") %>% 
    pivot_longer(-RowWell, names_to = "ColWell", values_to = "concentration") 
  
  
  pm_dilutions <- read.xlsx(pm_location, rows = c(32:40), cols = c(1:13), 
                            colNames = T, rowNames = T) %>% 
    rownames_to_column("RowWell") %>% 
    pivot_longer(-RowWell, names_to = "ColWell", values_to = "dilution") 
  
  pm_timepoint <- read.xlsx(pm_location, rows = c(44:52), cols = c(1:13), 
                            colNames = T, rowNames = T) %>% 
    rownames_to_column("RowWell") %>% 
    pivot_longer(-RowWell, names_to = "ColWell", values_to = "timepoint") 
  
  pm_studynr <- read.xlsx(pm_location, rows = c(55:63), cols = c(1:13), 
                          colNames = T, rowNames = T) %>%   
    rownames_to_column("RowWell") %>% 
    pivot_longer(-RowWell, names_to = "ColWell", values_to = "studynr", values_transform = list(studynr=as.character)) 
  
  #Merge all data in one dataframe
  pct_results <- bind_cols(pr_results, 
                           sample_id = pm_sample_ids$sample_id, 
                           studynr = pm_studynr$studynr,
                           concentration = pm_concentrations$concentration, 
                           dilution = pm_dilutions$dilution, 
                           timepoint = pm_timepoint$timepoint)
  
  ## Calculate mean value of duplicate measurements
  #mutate a colum to create (mean)RLU column, then for all samples: If another sample with the same sample_id exists, 
  #add their mean RFU to new column. Otherwise, add their original RFU to the new column
  pct_meanRLU <- data.frame("RowWell" = NA, "ColWell" = NA, "RLU" = NA, "sample_id" = NA, "studynr" = NA, "concentration" = NA, "dilution" = NA, "timepoint" = NA, "meanRLU" = NA)
  for (i in 1:length(unique(pct_results$sample_id))){
    sampleID <- unique(pct_results$sample_id)[[i]]
    x <- filter(pct_results, sample_id == sampleID)
    if(nrow(x) > 1) {
      meanRLU <- mean(x$RLU)
    } else {
      meanRLU <- x$RLU
    }
    x <- mutate(x, "meanRLU" = meanRLU)
    pct_meanRLU <- bind_rows(pct_meanRLU, x)
  }
  pct_meanRLU <- filter(pct_meanRLU, !is.na(sample_id)) 
  
  #Subtract mean blank from all samples
  blank <- filter(pct_meanRLU, sample_id == "Blank")$meanRLU[[1]]
  pct_scaled <- pct_meanRLU %>% 
    mutate(scaledRLU = meanRLU - blank)
  
  #Fit standard curve
  fit <- lm(scaledRLU~concentration, data=pct_scaled)
  intercept <- fit$coefficients[[1]]
  slope <- fit$coefficients[[2]]
  
  
  #Plot standard curve and add to return list
  CLIAstandard_curve <- ggplot(data = pct_scaled, aes(x = concentration, y = scaledRLU))+
    geom_point(na.rm = TRUE)+
    geom_abline(aes(intercept = intercept, slope = slope))+
    theme_bw()
  return_list[[1]] <- CLIAstandard_curve
  
  #Calculate concentrations using standard curve: RLU = intercept + [conc]*slope 
  # concentration = (RLU - intercept) / slope
  pct_conc <- pct_scaled %>% 
    mutate(PCT_sample_conc = (scaledRLU - intercept)/slope) %>% 
    mutate(PCT_concentration = PCT_sample_conc * dilution) %>% 
    mutate(rounded = round(PCT_concentration, digits = 2)) 
  
  # # put in original plate map format for readability
  # pct_concentrations <- pct_conc %>% 
  #   select(RowWell, ColWell, rounded) %>% 
  #   pivot_wider(names_from = ColWell, values_from = rounded) %>% 
  #   column_to_rownames(var = "RowWell") 
  
  # Select only useful columns from concentrations data: sample_id, studynr and rounded PCT concentration 
  pct_reduced <- pct_conc %>% 
    filter(is.na(concentration)) %>% 
    select(sample_id, studynr, rounded) %>% 
    distinct() %>% 
    filter(sample_id != c("Blank", "QC"))
  colnames(pct_reduced) <- c("sample.id.S", "subject.id", "pct")
  
  ## Replace concentration by "<6,86" if concentration is below < 6,86 and by ">5000" if concentration is over >5000
  pct_reduced <- pct_reduced %>% 
    mutate(pct_string = as.character(pct))
  for (i in 1:nrow(pct_reduced)){
    if (pct_reduced$pct[i] < 6.86) {pct_reduced$pct_string[i] <- "<6.86"} 
    if (pct_reduced$pct[i] > 5000) {pct_reduced$pct_string[i] <- ">5000"} 
  }
  
  return_list[[2]] <- pct_reduced
  
  return(return_list)
}

#Add CRP, PCT and creatinine to the data
AddCRPAndPCT <- function(dat) {
  # Load additional clinical patient data over time: temperature, CRP, leukocyte count, creatinine, information about antibiotic treatment 
  infl.markers.longitudinal <- read.spss("data/raw/aanvullende data ilona 1762021.sav", 
                                         use.value.labels = TRUE, to.data.frame = TRUE)
  
  # Make dataframe with CRP data over time
  crp <- select(infl.markers.longitudinal, "Studienr", "CRP_0", "CRP_dg1", "CRP_dg2", "CRP_dg4", "CRP_30") 
  names(crp) <- c("subject.id", "0", "1", "2", "4", "30")
  crp_long <- pivot_longer(crp, names_to = "day", values_to = "crp", -"subject.id", names_transform = list(day = as.integer))
  crp_long$subject.id <- as.character(crp_long$subject.id)
  
  # Make dataframe with creatinine data over time
  crea <- select(infl.markers.longitudinal, "Studienr", "kreat_dag0", "Kreat_1", "Kreat_2", "Kreat_4", "Kreat_30") 
  names(crea) <- c("subject.id", "0", "1", "2", "4", "30")
  crea_long <- pivot_longer(crea, names_to = "day", values_to = "crea", -"subject.id", names_transform = list(day = as.integer))
  crea_long$subject.id <- as.character(crea_long$subject.id)
  
  # Make dataframe with PCT data over time
  # Calculate PCT levels for all batches using function PCTcalculator
  # Batch1
  pr_results_B1 <- read.csv2("data/pct/CLIA_batch1_21-01-2022_gain3500.csv", skip = 8, sep = ":")
  pm_location_B1 <- "data/pct/plate_design_batch1.xlsx"
  pct_B1 <- PCTcalculator(pr_results_B1, pm_location_B1)
  # Batch2
  pr_results_B2 <- read.csv2("data/pct/CLIA_batch2_24-01-2022_gain3500.csv", skip = 8, sep = ":")
  pm_location_B2 <- "data/pct/plate_design_batch2.xlsx"
  pct_B2 <- PCTcalculator(pr_results_B2, pm_location_B2)
  # Batch3
  pr_results_B3 <- read.csv2("data/pct/CLIA_batch3_25-01-2022_gain3500.csv", skip = 8, sep = ":")
  pm_location_B3 <- "data/pct/plate_design_batch3.xlsx"
  pct_B3 <- PCTcalculator(pr_results_B3, pm_location_B3)
  # Merge results from the three batches in one dataframe:
  pct_df <- bind_rows(pct_B1[[2]], pct_B2[[2]], pct_B3[[2]])
  # MANUALLY adjust name of 1 sample because it is incomplete
  pct_df$sample.id.S[103] <- "303, S341185, 17-09-2010"
  # Add day information
  dayinfo <- select(dat, sample.id.S, subject.id, day) %>% 
    mutate(subject.id = as.character(subject.id))
  pct_long <- left_join(pct_df, dayinfo, by = c("sample.id.S", "subject.id"))
  
  AB_switch <- infl.markers.longitudinal %>% 
    mutate(subject.id = as.character(Studienr)) %>% select(subject.id, AB_switch) %>% 
    mutate(AB_switch = recode(AB_switch, ja = "Yes", nee = "No"))
  
  ## Add CRP and PCT to reduced data frame
  data.full <- dat %>% 
    left_join(crp_long, by = c("subject.id", "day")) %>%
    select(-CRP) %>% # remove incomplete CRP (day 0 only)
    left_join(pct_long, by = c("subject.id", "day", "sample.id.S")) %>% 
    left_join(AB_switch, by = "subject.id") %>% 
    left_join(crea_long, by =  c("subject.id", "day"))
  
  return(data.full)
}


#calculates different sums and ratios
AddRatiosAndSums <- function(data.clean) {
  ratios <- with(data.clean, 
      data.frame(sample.id = sample.id, pathogen = pathogen) %>% 
      #Add Sums
      mutate(BCAA_sum = L.Isoleucine + L.Leucine + L.Valine,
             TCA_cycle_sum = OA03_._Citric_acid + OA07_._Lactic_acid + OA08_._Malic_acid + OA12_._Fumaric_acid,
             urea_cycle_sum = Citrulline + L.Arginine + Ornithine + OA12_._Fumaric_acid,
             lc_Carnitines_sum = Myristoilcarnitine + Hexadecenoylcarntine + 
             Palmitoylcarnitine + Stearoylcarnitine + Dodecenoylcarnitine + Tetradecenoylcarnitine + 
               Linoleylcarnitine + Oleylcarnitine + Tetradecadienylcarntine,
             mc_Carnitines_sum = Hexanoylcarnitine + Octanoylcarnitine + Octenoylcarnitine + 
               Decanoylcarnitine + Lauroylcarnitine + Nonaylcarnitine + Pimeylcarnitine + Decenoylcarnitine,
             sc_Carnitines_sum = Acetylcarnitine + Propionylcarnitine + Isobutyrylcarnitine + 
               Butyrylcarnitine + Tiglylcarnitine + X2.methylbutyroylcarnitine + Isovalerylcarnitine,
             Cer_sum = Cer.d18.1.22.1. + Cer.d18.1.24.1. + Cer.d18.1.24.0. + Cer.d18.1.16.0. + 
               Cer.d18.1.23.0. + Cer.d18.0.24.0.,
             SM_sum = SM.d18.1.14.0. + SM.d18.1.15.0. + SM.d18.1.16.1. + SM.d18.1.16.0. + SM.d18.1.17.0. + 
               SM.d18.1.18.2. + SM.d18.1.18.1. + SM.d18.1.18.0. + SM.d18.1.20.1. + SM.d18.1.20.0. +
               SM.d18.1.21.0. + SM.d18.1.22.1. + SM.d18.1.22.0. + SM.d18.1.23.1. + SM.d18.1.23.0. + SM.d18.1.24.2. +
               SM.d18.1.24.1. + SM.d18.1.24.0. + SM.d18.1.25.1. + SM.d18.1.25.0.,
             LPC_sum = LPC.14.0. + LPC.16.0. + LPC.16.1. + LPC.18.0. + LPC.18.1. + LPC.18.2. + 
               LPC.18.3. + LPC.20.4. + LPC.20.5. + LPC.22.6. + LPC.O.16.1. + LPC.O.18.1.,
             PC_sum = PC.32.2. + PC.32.1. + PC.32.0. + PC.O.34.3. + PC.O.34.2. + PC.O.34.1. + PC.34.4. + 
               PC.34.3. + PC.34.2. + PC.34.1. + PC.O.36.6. + PC.O.36.5. + PC.O.36.4. + PC.O.36.3. + PC.O.36.2. + 
               PC.36.6. + PC.36.5. + PC.36.4. + PC.36.3. + PC.36.2. + PC.36.1. + PC.O.38.7. + PC.O.38.6. + 
               PC.O.38.5. + PC.O.38.4. + PC.38.7. + PC.38.6. + PC.38.5. + PC.38.4. + PC.38.3. + PC.38.2. + 
               PC.O.40.6. + PC.40.8. + PC.40.7. + PC.40.6. + PC.40.4. + PC.O.42.6. + PC.O.44.5. + PC.40.5.) %>% 
      #Add ratios
      mutate(HT5_Trp_ratio = Serotonine / L.Tryptophan,
             ADMA_Arg_ratio = ADMA / L.Arginine,
             SDMA_Arg_ratio = SDMA / L.Arginine,
             Carnitine_sum_lc_Carnitines_ratio = Carnitine / lc_Carnitines_sum,
             Carnitine_sum_mc_Carnitines_ratio = Carnitine / mc_Carnitines_sum, 
             Carnitine_sum_sc_Carnitines_ratio = Carnitine / sc_Carnitines_sum, 
             DCA_CA_ratio = DCA / CA,
             FA_14.1_14.0 = FA..14.1. / FA..14.0.,
             FA_16.1_16.0 = FA..16.1. / FA..16.0.,
             Gln_Glu = L.Glutamine / L.Glutamic.acid,
             Kyn_Trp = L.Kynurenine / L.Tryptophan,
             sum_BCAA_sum_Phe_Tyr_ratio =  BCAA_sum / (L.Phenylalanine + L.Tyrosine),
             sum_CER_sum_SM_ratio = Cer_sum / SM_sum,
             sum_LPC_sum_PC_ratio = LPC_sum / PC_sum))
  
  data.ratios <- data.clean %>% left_join(ratios, by = c("sample.id", "pathogen"))
  
  attr(data.ratios, "ratiorange") <- setdiff(colnames(ratios), c("sample.id", "pathogen"))
  
  return(data.ratios)
}

