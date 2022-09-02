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

#Save loadings
write.csv(pca_all_times[["rotation"]], file = "results/pca_loadings.csv", row.names = T)

# Visualize results
pca_data <- data.frame(data.pretreated[rownames(pca_metabolites), c("day", "subject.id", "hospitalization.time", "curb", "age", "sex")], 
                       pca_all_times$x[, 1:4])

plot_pca_per_timepoint <- pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = hospitalization.time)) +
  geom_point() +
  scale_color_lei(palette = "gradient", discrete = FALSE) +
  facet_grid(~ day) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom")

plot_pca_per_timepoint_begin <- pca_data %>% 
  mutate(day = paste("Day", day)) %>% 
  mutate(day = factor(day, levels = paste("Day", c(0, 1, 2, 4, 30)))) %>% 
  ggplot(aes(x = PC1, y = PC2, color = subject.id)) +
  geom_point(show.legend = F) +
  scale_color_lei(palette = "nine", discrete = TRUE) +
  facet_grid(~ day) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom")


plot_pca_per_timepoint_curb <- pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(curb))) +
  geom_point() +
  scale_color_lei(palette = "mixed", discrete = TRUE) +
  facet_grid(~ day) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom")

plot_PC_time_curve <- function(PC = "PC1", color = "subject.id") {
  pca_data %>% 
    ggplot(aes_string(x = "day", color = color, y = PC)) +
    geom_point(show.legend = F) +
    geom_line(aes(group = subject.id), show.legend = F) +
    scale_x_time_squish()  +
    theme_bw()
}

pca_plot <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot() +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1), show.legend = F, color = lei_colors["orange"]) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2), show.legend = F, color = lei_colors["orange"]) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4), show.legend = F, color = lei_colors["orange"]) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30), show.legend = F, color = lei_colors["orange"]) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F, color = lei_colors["black"]) +
  facet_wrap(~ subject.id) +
  labs(x = "PC1", y = "PC2") +
  coord_equal() +
  theme_bw()

pca_plot_curb <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot(aes(group = subject.id, color = as.factor(curb))) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F) +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1)) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2)) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4)) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30)) +
  theme_bw()

pca_plot_hosp_time <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot(aes(group = subject.id, color = hospitalization.time)) +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1)) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2)) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4)) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30)) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F) +
  scale_color_viridis_c() +
  theme_bw()

pca_plot_facetted <- pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = subject.id)) +
  geom_point(show.legend = F) +
  facet_grid(~ day) +
  theme_bw()


pca_plot_age <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot(aes(group = subject.id, color = age)) +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1)) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2)) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4)) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30)) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F) +
  scale_color_viridis_c() +
  theme_bw()

pca_plot_sex <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot(aes(group = subject.id, color = sex)) +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1)) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2)) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4)) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30)) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F) +
  theme_bw()


#Gather figures
figure6a <- plot_pca_per_timepoint_begin
figure6b <- plot_PC_time_curve("PC1", color = "hospitalization.time") +
  scale_color_lei(palette = "gradient", discrete = FALSE)
figure6c <- plot_PC_time_curve("PC2", color = "hospitalization.time") +
  scale_color_lei(palette = "gradient", discrete = FALSE)
figureS6 <- pca_plot

save(figure6a, figure6b, figure6c, figureS6, file = "manuscript/figures/plots_PCA.Rdata")
