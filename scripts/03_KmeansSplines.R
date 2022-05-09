#start with dat from 01_DataAnalysis
library(tidyverse)

spline_data <- as.data.frame(matrix(ncol = 5, nrow = length(metrange)))
rownames(spline_data) <- metrange

for (met in metrange) {
  spline_met <- smooth.spline(x = dat$day, y = dat[, met])
  spline_data[met, ] <- spline_met$y
}

k <- 5
kmeans_result <- kmeans(spline_data[, 1:5], k)

spline_data$cluster <- kmeans_result$cluster

colnames(spline_data)[1:5] <- sort(unique(dat$day))

long_cluster <- spline_data %>% 
  rownames_to_column("metabolite") %>% 
  pivot_longer(-c(metabolite, cluster), names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))

centroids <- as.data.frame(kmeans_result$centers) %>% 
  rownames_to_column("cluster") %>% 
  pivot_longer(-cluster, names_to = "day", values_to = "centroid") %>% 
  mutate(day = as.numeric(day))

long_cluster %>% 
  ggplot(aes(x = day, color = as.factor(cluster))) +
  geom_point(aes(y = spline)) +
  geom_line(aes(group = metabolite, y = spline), alpha = 0.1) +
  geom_line(data = centroids, aes(y = centroid), size = 2) +
  
  #geom_smooth(aes(group = as.factor(cluster)), se = F) +
  theme_bw()


#check number of metabolites per cluster
table(spline_data$cluster)


