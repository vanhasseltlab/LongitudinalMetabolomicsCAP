#04_PCT_measurement
library(tidyverse)
library(openxlsx)
library(ggplot2)

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
  colnames(pct_reduced) <- c("sample_id", "studynr", "pct")
  
  return_list[[2]] <- pct_reduced
    
  return(return_list)
  }

# Run function for all batches
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

## - Replace concentration by "<6,86" if concentration is below < 6,86
## - Replace concentration by ">5000" if concentration is over >5000
pct_LDs <- pct_df %>% 
  mutate(pct_string = as.character(pct))
for (i in 1:nrow(pct_LDs)){
  if (pct_LDs$pct[i] < 6.86) {pct_LDs$pct_string[i] <- "<6.86"} 
  if (pct_LDs$pct[i] > 5000) {pct_LDs$pct_string[i] <- ">5000"} 
}


## TODO: 
## - Make figure with lines through points per patient, x time, y concentration. 

#Plot samples over time
sampleset <- filter(pct_conc, !is.na(timepoint))
standards <- filter(pct_conc, !is.na(concentration)) %>% 
  mutate(timepoint = 0)
ggplot(data = sampleset)+
  geom_point(aes(x = timepoint, y = rounded, color = as.factor(studynr)), size = 2, alpha = 0.6)+
  geom_point(data = standards, aes(x = timepoint, y = concentration), color = "black")+
  # scale_y_log10()+
  theme_bw()

