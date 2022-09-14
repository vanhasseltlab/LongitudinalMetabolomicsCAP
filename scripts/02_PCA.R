# Exploratory analysis: visualize trends through dimension reduction (PCA) 

# Load libraries
library(tidyverse)

# Load functions
source("functions/LeidenColoring.R")

# Load data from 01_DataPreparation (data.pretreated)
load("data/01_data_clean.Rdata")

# Select only metabolites
metrange <- attr(data.pretreated, "metrange")
metrange <- metrange[colSums(is.na(data.pretreated[, metrange])) < 2]

pca_metabolites <- data.pretreated[!is.na(rowSums(data.pretreated[, metrange])), metrange]

# PCA
set.seed(123)
pca_all_times <- prcomp(pca_metabolites)
percentage_explained <-  round(summary(pca_all_times)$importance[2, 1:2]*100, 1)

#Save loadings
write.csv(pca_all_times[["rotation"]], file = "results/pca_loadings.csv", row.names = T)

# Visualize results
pca_data <- data.frame(data.pretreated[rownames(pca_metabolites), c("day", "subject.id", "hospitalization.time", "curb")], 
                       pca_all_times$x[, 1:4])

plot_pca_per_timepoint <- pca_data %>% 
  mutate(day = paste("Day", day)) %>% 
  mutate(day = factor(day, levels = paste("Day", c(0, 1, 2, 4, 30)))) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(show.legend = F, color = "black", alpha = 0.7, stroke = 0, size = 2) +
  facet_grid(~ day) +
  labs(x = paste0("PC1 (", percentage_explained[1], "%)"), 
       y = paste0("PC2 (", percentage_explained[2], "%)")) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill = "white"))

plot_pca_individuals <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot() +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1), show.legend = F, color = lei_colors["soft orange"]) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2), show.legend = F, color = lei_colors["soft orange"]) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4), show.legend = F, color = lei_colors["soft orange"]) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30), show.legend = F, color = lei_colors["soft orange"]) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F, color = lei_colors["black"]) +
  facet_wrap(~ subject.id) +
  labs(x = "PC1", y = "PC2") +
  coord_equal() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))


#Gather figures
figure1 <- plot_pca_per_timepoint
figureS1 <- plot_pca_individuals

save(figure1, figureS1, file = "manuscript/figures/plots_PCA.Rdata")
