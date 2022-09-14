### LongitudinalMetabolomicsCAP

The scripts and data required for replication of the study 'Metabolic profiles over time for treatment response monitoring in community-acquired pneumonia'

#### Data

Anonymized data was added in the folder `data/`. The file `data/01_data_clean.Rdata` contains the metabolomics measurements, inflammatory markers and CURB on different time points for each of 25 included patients. This data can be used to run the scripts from `02_PCA.R` onwards.

#### Scripts (`scripts/`)

-   `00_DataAssembly.R`, combines patient data with metabolomics measurements from the different platforms to one large dataset 'data/00_data_raw.csv'

-   `01_DataPreparation.R`, prepares data through filtering, selecting, scaling and calculating new values using the functions from `AllFunctions.R`. Creates `data/01_data_clean.Rdata`

-   `02_PCA.R`, use principal component analysis and visualize PC's

-   `03_KmeansSplines.R`, capture each metabolite over all patients in a spline and cluster similar metabolite profiles using Kmeans clustering

-   `04_Correlations.R`, calculate correlations between the (change in) metabolites and PCT, CRP, length of stay and CURB

-   `05_FiguresManuscript.R`, loads all ggplot objects generated in the previous scripts to the environment and creates png files for each figure for the manuscript.

#### Helper function scripts (`functions/`)

-   `AllFunctions.R`, contains all functions used for data preparation

-   `LeidenColoring.R`, contains color used in figures
