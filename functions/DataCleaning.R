DataCleaning <- function(dat, metrange, nonmetrange) {
  # This function removes samples (rows) with no metabolomics measurements in at 
  # least one platform (>10 missing metabolites) and then removes metabolites 
  # (columns) with missing values.
  
  # Subset metabolites from data
  dat.subset <- dat[, metrange]
  # Calculate the numer of NA's per row
  narows <- apply(dat.subset, 1, function(X){sum(is.na(X))})
  # Subset metabolomics data keeping only rows with les then 10 missing metabolites
  dat.2 <- dat[narows < 10, ]
  
  # Subset metabolites from partly cleaned dataset
  dat.subset.2 <- dat.2[, metrange]
  # Calculate the numer of NA's per column
  nacols <- apply(dat.subset.2, 2, function(X){sum(is.na(X))})
  # keep column names of metabolites without missing metabolites
  keep <- colnames(dat.subset.2[, nacols==0])
  
  # Combine cleaned metabolomicsdat with non-metabolomicsdat 
  #  and return clean dat set.
  dat.clean <- cbind(dat.2[, keep], dat.2[, nonmetrange])
  return(dat.clean)
}