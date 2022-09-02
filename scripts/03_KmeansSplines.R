#Fit splines and Kmeans to cluster metabolite curves

# Load libraries
library(tidyverse)
library(readxl)

# Source functions
source("functions/LeidenColoring.R")

# Load data from 01_DataPreparation (data.full.clean, data.pretreated)
load("data/01_data_clean.Rdata")

### Calculate splines for each metabolite, over all patients
metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))

# Calculate splines for each metabolite
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
  geom_line(data = centroids, aes(y = centroid), size = 1.3) +
  scale_color_lei(palette = "five") +
  labs(color = "Cluster", x = "Time (days)", y = "Splines over patients\n(scaled metabolite values)") +
  theme_bw()

#check number of metabolites per cluster
table(spline_data$cluster)

###### Metabolic classes and clusters ####
names_and_class <- read.csv("data/metabolite_names_and_classes.csv", stringsAsFactors = F)

metabolite_class <- long_cluster %>%
  select(metabolite, cluster) %>% 
  distinct()  %>%  left_join(names_and_class)

write.csv(metabolite_class, file = "results/kmeans_clusters_metabolites.csv", row.names = F)

#### Visualize metabolite classes per cluster #####
plot_metabolite_class <- metabolite_class %>% 
  filter(!is.na(class)) %>% 
  left_join(metabolite_class) %>% 
  ggplot(aes(x = class, fill = as.factor(cluster))) +
  geom_bar() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_lei(palette = "five") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Metabolite class", y = "Count", fill = "Cluster") +
  theme_bw()

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
      filter(!is.na(value)) %>% 
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

plot_patient_clusters <- plot_data %>% 
  ggplot(aes(x = day, y = value, group = cluster, color = as.factor(cluster))) +
  geom_hline(yintercept = 0, color = "grey65") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, color = "grey65") +
  geom_line() +
  geom_point() +
  facet_wrap(~ subject.id, scales = "free_y") +
  scale_color_lei(palette = "five") +
  scale_x_time_squish() +
  labs(x = "Time (days)", y = "Cluster of scaled metabolite values", color = "Cluster") +
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
long_crp <- data.frame(day = spline_crp$x, spline = spline_crp$y) %>% 
  left_join(dat_crp_long %>% group_by(day) %>% dplyr::summarize(mean_crp = mean(crp_value, na.rm = T)))

#### Visualize CRP and clustered metabolite data ####
#metabolites_of_interest <- unique(c("PC.34.1.", grep("LPC.", long_cluster$metabolite, value = TRUE)))
#4:0, 16:0, 16:1, 18:3, en 20:4
metabolites_of_interest <- c("PC.34.1.", "LPC.14.0.", "LPC.16.0.", "LPC.16.1.", "LPC.18.3.", "LPC.20.4.")
metabolite_names <- read.xlsx("data/Metabolite_names_V2_incl_ratios.xlsx") %>% 
  dplyr::rename(metabolite = Detected_metabolite_name_in_R, metabolite_name = Metabolite_name_V2) %>% 
  select(-Metabolite_name)

plot_data <- long_cluster %>% 
  filter(metabolite %in% metabolites_of_interest) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  bind_rows(long_crp %>%  mutate(metabolite = "CRP", cluster = "CRP") %>% select(metabolite, cluster, day, spline)) %>% 
  group_by(metabolite, cluster) %>%
  dplyr::mutate(y_label = spline[day == 30]) %>% 
  ungroup() %>% 
  dplyr::mutate(CRP = ifelse(cluster == "CRP", "CRP", "")) %>% 
  left_join(metabolite_names) %>% 
  mutate(metabolite = ifelse(is.na(metabolite_name), metabolite, metabolite_name))

lpc_plot <- plot_data %>% 
  ggplot(aes(x = day, color = cluster)) +
  geom_line(aes(group = metabolite, y = spline, linetype = cluster), alpha = 0.7) +
  geom_text(aes(label = metabolite, x = 30, y = y_label), hjust = -0.02, show.legend = F, size = 2) +
  scale_color_lei(palette = "six") +
  scale_linetype_manual(values = c(rep(1, 2), 2)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Time (days)", y = "Scaled metabolite and CRP levels", color = "Cluster", linetype = "Cluster") +
  theme_bw()

#### Instead of cluster use metabolite class for splines ####

data.class.long <- data.clusters.long %>% 
  left_join(names_and_class) %>% filter(!is.na(class)) %>% 
  filter(!is.na(value))


nr_classes <- length(unique(data.class.long$class))
spline_data_class <- as.data.frame(matrix(ncol = 5, nrow = nr_classes))
rownames(spline_data_class) <- unique(data.class.long$class)
for (met_class in unique(data.class.long$class)) {
  dat <- data.class.long %>% 
    filter(class %in% met_class)
  spline_met <- smooth.spline(x = dat$day, y = dat$value)
  spline_data_class[met_class, ] <- spline_met$y
}


colnames(spline_data_class)[1:5] <- sort(unique(dat$day))
long_class <- spline_data_class %>% 
  rownames_to_column("class") %>% 
  pivot_longer(-class, names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))


### Visualization ####
table_classes <- sort(table(metabolite_class$class))
small_classes <- names(table_classes[table_classes < 3])

cluster_plot_class_facets <- long_cluster %>% 
  left_join(metabolite_class) %>% 
  filter(!is.na(class)) %>% 
  filter(!class %in% small_classes) %>% 
  mutate(class = gsub("-", "- ", class)) %>% 
  mutate(class = gsub("phosphatidyl", "phosphatidyl ", class)) %>% 
  ggplot(aes(x = day)) +
  geom_line(aes(y = spline, group = metabolite, color = as.factor(class)), alpha = 0.5) +
  geom_line(data = long_class %>% filter(!class %in% small_classes) %>% 
              mutate(class = gsub("-", "- ", class)) %>% 
              mutate(class = gsub("phosphatidyl", "phosphatidyl ", class)), 
            aes(y = spline), color = "black", show.legend = F) +
  facet_wrap(~ class, labeller = labeller(class = label_wrap_gen(15))) +
  scale_x_time_squish() +
  theme_bw() +
  theme(legend.position = "none") 


#Gather figures
figure2 <- cluster_plot
figureS1 <- plot_metabolite_class
figureS2 <- cluster_plot_class_facets
figureS3 <- plot_patient_clusters
figure3b <- lpc_plot
tableS2 <- select(metabolite_class, -metabolite)


save(figure2, figureS1, figureS2, figureS3, figure3b, tableS2, file = "manuscript/figures/plots_Kmeans.Rdata")
