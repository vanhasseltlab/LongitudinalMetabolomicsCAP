## This script contains all functions nessecary to reproduce the data from the Longitudinal Metabolomics project

#### Remove redundant variables from raw data using reduction function ----
ReduceData <- function(data.raw){
  metabolomicsdata <- data.raw[, 6:374]
  patientdata <- select(data.raw, 
                        "sample.id"= NMC.name, # Metabolomics sample id
                        "sample.id.S" = SampleID, # Hospital sample id
                        "subject.id" = Studynr, # Study patient id
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
  metadata <- c("sample.id", "sample.id.S", "subject.id", "day", "pathogen", "pathogen.group", "age", "sex", "psi.score", "nursing.home.resident",
                "renal.disease", "liver.disease", "congestive.heart.failure", "cns.disease", "malignancy", "altered.mental.status",
                "respiratory.rate", "systolic.blood.pressure", "temperature", "pulse", "pH" , "BUN" , "sodium", "glucose", 
                "hematocrit", "partial.pressure.of.oxygen", "pleural.effusion.on.x.ray", "race", "duration.of.symptoms.before.admission",
                "antibiotic.treatment.before.admission", "corticosteroid.use.before.admission", "COPD", "diabetes", 
                "oxygen.saturation", "supplemental.oxygen.required", "antibiotics.started.in.hospital", "CRP", "leukocyte.count",
                "IL.6", "IL.10")
  
  metrange <- setdiff(names(data.reduced), metadata)
  
  # Define metabolite data as numeric data for analysis
  data.reduced[,metrange] <- apply(data.reduced[,metrange], 2, function(x) {
    as.numeric(as.character(x))
  })
  
  return(data.reduced)
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





#### 
#### Custom color palette function ----
# Create vector with 8 Leiden university / LACDR colors
lei_colors <- c(
  `blue`       = "#001158",
  `orange`     = "#FF9933", 
  `red`        = "#be1908",
  `lightgreen` = "#aaad00",
  `brightgreen`= "#06bd09",
  `darkgreen`  = "#2c712d",
  `turquoise`  = "#34a3a9",
  `lightblue`  = "#5cb1eb",
  `brightblue` = "#0536fa",
  `violet`     = "#b02079" )

## Function to extract the hex codes from this vector by name
#' Function to extract lei colors as hex codes
#' @param ... Character names of lei_colors 
lei_cols <- function (...){
  cols <- c(...)
  
  if (is.null(cols))
    return (lei_colors)
  
  lei_colors[cols]
}

lei_palettes <- list(
  `main`  = lei_cols("blue", "orange"),
  `three` = lei_cols("blue", "orange", "darkgreen"),
  `cool`  = lei_cols("blue", "lightblue", "turquoise", "lightgreen", "darkgreen"),
  `hot`   = lei_cols("violet", "red", "orange"),
  `mixed` = lei_cols("blue", "lightblue", "turquoise", "darkgreen", "lightgreen", "orange", "red", "violet"),
  `two`   = lei_cols("red", "violet"), 
  `five`  = lei_cols("blue", "lightblue", "orange", "red", "darkgreen"),
  `nine`  = lei_cols("lightblue", "violet", "brightgreen",  "brightblue", "red", "lightgreen", "blue", "orange", "darkgreen"))

#' Return function to interpolate a lei color palette
#' @param palette Character name of palette in lei_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette(), such as an alpha value
lei_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- lei_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

# Now return a function for any palette, for example `cool`:
lei_pal("cool")
# The returned function will interpolate the palette colors for a certain number of levels, making it possible to create shades between our original colors. To demonstrate, we can interpolate the "cool" palette to a length of 10:
lei_pal("cool")(10)
# This is what we need to create custom ggplot2 scales

## Create custom color and fill scales for ggplot2 by creating one function for color and one for fill. 
#' Color scale constructor for lei colors
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for lei colors
#'
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}





#### 