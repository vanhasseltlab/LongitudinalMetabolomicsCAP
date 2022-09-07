### Scripts for 'Metabolic profiles over time for treatment response monitoring in community-acquired pneumonia'

##### Scripts

-   `00_DataAssembly.R`, combines patient data with metabolomics measurements from the different platforms to one large dataset 'data/00_data_raw.csv'

-   01_DataPreparation.R, prepares data through filtering, selecting, scaling and calculating new values using the functions from AllFunctions.R 

-   02_PCA.R, use principal component analysis and visualize PC's

-   03_KmeansSplines.R, capture each metabolite over all patients in a spline and cluster similar metabolite profiles using Kmeans clustering

-   04_Correlations.R, calculate correlations between the (change in) metabolites and PCT, CRP, length of stay and CURB

-   GeneralScript.R, can be used to source each script. (in progress)