#Fit splines and Kmeans to cluster metabolite curves

# Load libraries
library(tidyverse)
library(readxl)
# Source functions
source("functions/AllFunctions.R")
source("functions/LeidenColoring.R")

#Data preparation
dat.raw <- read.csv2("data/00_data_raw.csv")
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
#Remove metabolites with to many NAs
data.clean <- DataCleaning(data.reduced, attr(data.reduced, "metrange"))

#Add longitudinal markers: CRP and PCT and ratios and sums
data.pc <- AddCRPAndPCT(data.clean)
data.rs <- AddRatiosAndSums(data.pc)

#Scale metabolite values
data.pretreated <- DataPretreatment(data.rs, metrange = c(attr(data.rs, "metrange"), attr(data.rs, "ratiorange")))

### Calculate splines for each metabolite, over all patients
metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))

spline_data <- as.data.frame(matrix(ncol = 5, nrow = length(metrange)))
rownames(spline_data) <- metrange
for (met in metrange) {
  dat <- data.pretreated[, c(met, "day")] %>% 
    filter(!is.na(!!as.name(met)))
  spline_met <- smooth.spline(x = dat$day, y = dat[, met])
  spline_data[met, ] <- spline_met$y
}

### Use Kmeans to cluster similar metabolite profiles
set.seed(73485183)
k <- 5
kmeans_result <- kmeans(spline_data[, 1:5], k)

# Add results to data
spline_data$cluster <- kmeans_result$cluster
colnames(spline_data)[1:5] <- sort(unique(dat$day))
long_cluster <- spline_data %>% 
  rownames_to_column("metabolite") %>% 
  pivot_longer(-c(metabolite, cluster), names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))

colnames(kmeans_result$centers)[1:5] <- sort(unique(dat$day))
centroids <- as.data.frame(kmeans_result$centers) %>% 
  rownames_to_column("cluster") %>% 
  pivot_longer(-cluster, names_to = "day", values_to = "centroid") %>% 
  mutate(day = as.numeric(day))

### Visualize clusters ####
cluster_plot <- long_cluster %>% 
  ggplot(aes(x = day, color = as.factor(cluster))) +
  geom_point(aes(y = spline)) +
  geom_line(aes(group = metabolite, y = spline), alpha = 0.1) +
  geom_line(data = centroids, aes(y = centroid), size = 2) +
  scale_color_lei(palette = "five") +
  theme_bw()

pdf("results/figures/cluster_plot_splines_metabolites.pdf", width = 6, height = 4)
print(cluster_plot)
dev.off()

#check number of metabolites per cluster
table(spline_data$cluster)
#####

###### Metabolic classes and clusters ####
# up to date metabolite names: read_excel("data/Metabolite_names_V2_incl_ratios.xlsx")
metabolite_class_raw <- read_excel("data/Biochemical_classes_metabolites.xlsx")
metabolite_names <- read_excel("data/Metabolite_names_V2_incl_ratios.xlsx")

metabolite_clusters <- long_cluster %>% 
  select(metabolite, cluster) %>% 
  distinct()

metabolite_class <- metabolite_class_raw %>%
  left_join(metabolite_names %>% mutate(Metabolite_name = Metabolite_name_V2) %>% 
              select(Metabolite_name, Detected_metabolite_name_in_R)) %>% 
  dplyr::rename(metabolite = Detected_metabolite_name_in_R) %>% 
  right_join(metabolite_clusters) %>% 
  mutate(Metabolite_class = ifelse(metabolite %in% attr(data.pretreated, "ratiorange"),
                                   "Sums and ratios", Metabolite_class))

#### Visualize metabolite classes per cluster #####
table_classes <- sort(table(metabolite_class$Metabolite_class))
small_classes <- names(table_classes[table_classes < 9])

cluster_plot_metabolite_class <- long_cluster %>% 
  left_join(metabolite_class) %>% 
  # mutate(Metabolite_class = ifelse(Metabolite_class %in% small_classes | is.na(Metabolite_class), "Other", Metabolite_class)) %>% 
  ggplot(aes(x = day)) +
  geom_line(data = centroids, aes(y = centroid), size = 1) +
  geom_point(aes(y = spline, color = Metabolite_class)) +
  geom_line(aes(group = metabolite, y = spline, color = Metabolite_class), alpha = 0.1) +
  facet_wrap(~ cluster) +
  theme_bw()

plot_metabolite_class <- metabolite_class %>% 
  ggplot(aes(x = Metabolite_class, fill = Metabolite_class)) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(~ cluster) +
  theme_bw()

plot_metabolite_class_alt <- metabolite_class %>% 
  left_join(metabolite_class) %>% 
  ggplot(aes(x = Metabolite_class, fill = as.factor(cluster))) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_lei(palette = "five") +
  theme_bw()
#####

#### Clustered metabolite profiles per patient ####
##Splines per patient (get the profiles!)
data.clusters.long <- data.pretreated %>% 
  select(c(all_of(metrange), day, subject.id)) %>% 
  pivot_longer(-c(day, subject.id), names_to = "metabolite") %>% 
  left_join(spline_data %>% rownames_to_column("metabolite") %>% 
              select(cluster, metabolite)) %>% 
  arrange(day)

patients <- unique(data.clusters.long$subject.id)
clusters <- unique(data.clusters.long$cluster)

spline_patient <- as.data.frame(matrix(ncol = 5, nrow = length(patients)*length(clusters)))
names(spline_patient) <- paste0("V", c("0", "1", "2", "4", "30"))
spline_patient$subject.id <- rep(patients, length(clusters))
spline_patient$cluster <- rep(clusters, each = length(patients))
rownames(spline_patient) <- paste0(spline_patient$cluster, spline_patient$subject.id)

for (patient in patients) {
  for (cl in sort(clusters)) {
    combination <-  paste0(cl, patient)
    dat <- data.clusters.long %>% 
      filter(cluster == cl & subject.id == patient) %>% 
      select(day, subject.id, value)
    if (length(unique(dat$day)) < 4) {
      mean_measurement <- dat %>% group_by(day) %>% 
        dplyr::summarize(mean_value = mean(value))
      
      spline_patient[combination, paste0("V", mean_measurement$day)] <- mean_measurement$mean_value
      next
    }
    spline_met <- smooth.spline(x = dat$day, y = dat$value)
    spline_patient[combination, 1:5] <- predict(spline_met, x = c(0, 1, 2, 4, 30))$y
  }
}

#### Visualize cluster per patient ####
plot_data <- spline_patient %>% 
  pivot_longer(-c(subject.id, cluster), names_prefix = "V", names_to = "day") %>% 
  mutate(day = as.numeric(day)) %>%
  left_join(data.pretreated %>% select(subject.id, curb, age, hospitalization.time) %>% distinct()) %>% 
  filter(!is.na(value))

plot_data %>% 
  ggplot(aes(x = day, y = value, group = cluster, color = as.factor(cluster))) +
  geom_hline(yintercept = 0, color = "grey65") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, color = "grey65") +
  geom_line() +
  geom_point() +
  facet_wrap(~ subject.id, scales = "free_y") +
  scale_color_lei(palette = "five") +
  theme_bw()

plot_data %>% 
  filter(subject.id != "78") %>% 
  #filter(day != 30) %>% 
  ggplot(aes(x = day, y = value, color = as.factor(curb))) +
  geom_line(aes(group = subject.id), alpha = 0.7) +
  geom_point() +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_bw()

plot_data %>% 
  filter(day != 30) %>% 
  ggplot(aes(x = day, y = value, color = log(hospitalization.time))) +
  geom_line(aes(group = subject.id), alpha = 0.7) +
  geom_point() +
  scale_color_lei(discrete = F) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_bw()

#### Splines over all metabolites for each patient ####
spline_one <- as.data.frame(matrix(ncol = 5, nrow = length(patients)))
names(spline_one) <- paste0("V", c("0", "1", "2", "4", "30"))
spline_one$subject.id <- patients
rownames(spline_one) <- spline_one$subject.id

for (patient in unique(data.clusters.long$subject.id)) {
  dat <- data.clusters.long %>% 
    filter(subject.id == patient) %>% 
    select(day, subject.id, value)
  
  if (length(unique(dat$day)) < 4) {
    mean_measurement <- dat %>% group_by(day) %>% 
      dplyr::summarize(mean_value = mean(value))
    
    spline_one[patient, paste0("V", mean_measurement$day)] <- mean_measurement$mean_value
    next
  }
  spline_met <- smooth.spline(x = dat$day, y = dat$value)
  spline_one[patient, 1:5] <- predict(spline_met, x = c(0, 1, 2, 4, 30))$y
}

plot_data_one <- spline_one %>% 
  pivot_longer(-c(subject.id), names_prefix = "V", names_to = "day") %>% 
  mutate(day = as.numeric(day)) %>%
  left_join(data.pretreated %>% select(subject.id, curb, age) %>% distinct()) %>% 
  filter(!is.na(value))

plot_data_one %>%  
  ggplot(aes(x = day, group = subject.id, y = value, color = as.factor(curb))) +
  geom_line(alpha = 0.7) +
  geom_point() +
  theme_bw()


#### Spline CRP ####
infl.markers.longitudinal <- read.spss("data/raw/aanvullende data ilona 1762021.sav", 
                                       use.value.labels = TRUE, to.data.frame = TRUE)
dat_crp <- select(infl.markers.longitudinal, Studienr, CRP_0, CRP_dg1, CRP_dg2, CRP_3, CRP_dg4, CRP_dg5, CRP_dg6, CRP_7, CRP_dg10, CRP_30)
dat_crp_long <- dat_crp %>% 
  pivot_longer(-c(Studienr), names_to = "day", values_to = "crp_value") %>% 
  na.omit() %>% #Remove NAs
  mutate(crp_value = log2(crp_value + 1)) %>% #Log transform
  mutate(crp_value = (crp_value - mean(crp_value)) / sd(crp_value)) %>% #Scale
  mutate(day = as.numeric(gsub("CRP_|CRP_dg", "", day)))
spline_crp <- smooth.spline(x = as.numeric(as.character(dat_crp_long$day)), y = dat_crp_long$crp_value)
long_crp <- data.frame(day = spline_crp$x, spline = spline_crp$y)

#### Visualize CRP and clustered metabolite data ####
# Separate plot for CRP
long_crp %>% 
  ggplot(aes(x = day, y = spline))+
  geom_point(color = "black")+
  geom_line(color = "black")+
  geom_point(data = dat_crp_long, aes(x = as.numeric(as.character(day)), y = crp_value, color = as.factor(Studienr)))+
  scale_color_lei(palette = "mixed")+
  labs(x = "Days after hospital admission", y = "Scaled CRP levels", color = "Patient number")+
  theme_bw()


# Plotting
long_cluster %>% 
  ggplot(aes(x = day, color = as.factor(cluster))) +
  geom_point(aes(y = spline), size = 0.5) +
  geom_line(aes(group = metabolite, y = spline), alpha = 0.1) +
  geom_line(data = centroids, aes(y = centroid), size = 1, alpha = 0.9) +
  scale_color_lei(palette = "five") +
  geom_line(data = long_crp, aes(x = day, y = spline), size = 0.75, alpha = 0.9, color = "black", linetype = "dashed")+
  labs(x = "Days after hospital admission", y = "Scaled metabolite levels", color = "Cluster", linetype = "Cluster")+
  theme_bw()




