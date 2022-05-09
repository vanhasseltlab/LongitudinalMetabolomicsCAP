#start with dat from 01_DataAnalysis
library(tidyverse)
dat = data.pretreated
dat <- logdata
spline_data <- as.data.frame(matrix(ncol = 5, nrow = length(metrange)))
rownames(spline_data) <- c(metrange)

for (met in metrange) {
  spline_met <- smooth.spline(x = dat$day, y = dat[, met])
  spline_data[met, ] <- spline_met$y
}

k <- 5
kmeans_result <- kmeans(spline_data[, 1:5], k)

spline_data$cluster <- kmeans_result$cluster

colnames(spline_data)[1:5] <- sort(unique(dat$day))

# Calculating centroids
long_cluster <- spline_data %>% 
  rownames_to_column("metabolite") %>% 
  pivot_longer(-c(metabolite, cluster), names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))

centroids <- as.data.frame(kmeans_result$centers) %>% 
  rownames_to_column("cluster")
colnames(centroids) <- c("cluster", sort(unique(dat$day)))
centroids <- centroids %>% 
  pivot_longer(-cluster, names_to = "day", values_to = "centroid") %>%
  mutate(day = as.numeric(day))

# Spline CRP
dat_crp <- select(infl.markers.longitudinal, Studienr, CRP_0, CRP_dg1, CRP_dg2, CRP_3, CRP_dg4, CRP_dg5, CRP_dg6, CRP_7, CRP_dg10, CRP_30)
dat_crp_long <- dat_crp %>% 
  pivot_longer(-c(Studienr), names_to = "day", values_to = "crp_value") %>% 
  na.omit() %>% #Remove NAs
  mutate(crp_value = log2(crp_value + 1)) #%>% #Log transform
  # mutate(crp_value = (crp_value - mean(crp_value)) / sd(crp_value)) #Scale

dat_crp_long$day <- as.factor(dat_crp_long$day)
levels(dat_crp_long$day) <- c("0", "3", "30", "7", "1", "10", "2", "4", "5", "6")

spline_data_crp <- as.data.frame(matrix(ncol = 10, nrow = 1))
colnames(spline_data_crp) <- c("0", "1", "2", "3", "4", "5", "6", "7", "10", "30")
rownames(spline_data_crp) <- "crp"

spline_crp <- smooth.spline(x = as.numeric(as.character(dat_crp_long$day)), y = dat_crp_long$crp_value)
spline_data_crp[1,] <- spline_crp$y

long_crp <- spline_data_crp %>% 
  rownames_to_column("crp") %>%
  pivot_longer(-c(crp), names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))

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
  scale_color_lei(palette = "five")+
  # geom_point(data = long_crp, aes(x = day, y = spline), size = 0.5, color = "black")+
  # geom_line(data = long_crp, aes(x = day, y = spline), size = 0.75, alpha = 0.9, color = "black", linetype = "dashed")+
  labs(x = "Days after hospital admission", y = "Scaled metabolite levels", color = "Cluster", linetype = "Cluster")+
  theme_bw()
ggsave("cluster_plot.png", width = 6, height = 4)

#plotting

#check number of metabolites per cluster
table(spline_data$cluste)

