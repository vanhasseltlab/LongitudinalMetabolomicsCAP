# Exploratory analysis: visualize trends through dimension reduction (PCA) 

# Load libraries
library(gridExtra)
library(tidyverse)
library(missMDA)

# Load functions
source("functions/LeidenColoring.R")

# Load data from 01_DataPreparation (data.pretreated)
load("data/01_data_clean.Rdata")
names_and_class <- read.csv("data/metabolite_names_and_classes.csv", stringsAsFactors = F)

# Select only metabolites
metrange <- attr(data.pretreated, "metrange")

#### deal with missingness through missMDA package
pca_metabolites <- imputePCA(data.pretreated[, metrange], method = "EM", ncp = 2)$completeObs

# PCA
pca_all_times <- prcomp(pca_metabolites)
percentage_explained <-  round(summary(pca_all_times)$importance[2, 1:2]*100, 1)

#Save loadings
write.csv(pca_all_times[["rotation"]], file = "results/pca_loadings.csv", row.names = T)

# Results
pca_data <- data.frame(data.pretreated[, c("day", "subject.id", "hospitalization.time", "curb")], 
                       pca_all_times$x[, 1:2]) %>% 
  mutate(PC1 = -PC1) #mirror PC1

loadings <- data.frame(metabolite = colnames(pca_metabolites), t(t(pca_all_times$rotation[, 1:2]) * pca_all_times$sdev[1:2]), stringsAsFactors = F) %>% 
  mutate(PC1 = -PC1) %>% #mirror PC1
  mutate(weights = sqrt(PC1^2 + PC2^2)) %>% 
  left_join(names_and_class) %>% 
  arrange(class, desc(weights)) %>% 
  mutate(var_explained_PC1 = PC1^2/sum(PC1^2),
         var_explained_PC2 = PC2^2/sum(PC2^2))

#variance explained by time component
lm_time <- lm(day ~ poly(PC1, 2, raw = TRUE) + poly(PC2, 2, raw = TRUE), data = pca_data)
variance_explained <- summary(lm_time)$adj.r.squared
lm_time

#variance explained by each class
class_var_explained <- loadings %>% group_by(class) %>% 
  summarize(var_explained_PC1 = sum(var_explained_PC1),
            var_explained_PC2 = sum(var_explained_PC2))

class_var_explained %>% 
  pivot_longer(-class, names_to = "PC", values_to = "expl_var", names_prefix = "var_explained_") %>% 
  ggplot(aes(y = class, x = expl_var)) +
  geom_point() +
  facet_grid(~ PC) +
  theme_bw()

# cholesteryl esters, LPC's, sphingomyelins, diacylglycerols, and triglycerides 
class_var_explained %>% 
  filter(class %in% c("Cholesteryl esters", "Lysophospholipids", "Sphingomyelin", "Diacylglycerols", "Triglycerides")) %>% 
  select(-class) %>% colSums %>% mean()*100

##more plots
pca_data %>% 
  ggplot(aes(x = PC1, y = rank(day))) +
  geom_point() +
  theme_bw()

pca_data %>% 
  ggplot(aes(x = PC2, y = rank(day))) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw()

pca_data %>% 
  ggplot(aes(x = 0.98376*PC1 + 0.05026*PC1^2, y = day)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw()


lm_time <- lm(rank(day) ~ PC1 + PC2, data = pca_data)
variance_explained <- summary(lm_time)$adj.r.squared
lm_time

# Visualize results
plot_pca_per_timepoint <- pca_data %>% 
  mutate(day = paste("Day", day)) %>% 
  mutate(day = factor(day, levels = paste("Day", c(0, 1, 2, 4, 30)))) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(show.legend = F, alpha = 0.7, stroke = 0, size = 2) +
  facet_grid(~ day) +
  labs(x = paste0("PC1 (", percentage_explained[1], "%)"), 
       y = paste0("PC2 (", percentage_explained[2], "%)")) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_rect(fill = "white"))

plot_pca_individuals <- pca_data %>% 
  left_join(data.pretreated %>% select(subject.id, day, hospitalization.time, curb)) %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2)) %>% 
  ggplot(aes(color = "#8592BC")) +
  geom_hline(yintercept = 0, color = "grey65") +
  geom_vline(xintercept = 0, color = "grey65") +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1), show.legend = F) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2), show.legend = F) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4), show.legend = F) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30), show.legend = F) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F, color = lei_colors["black"]) +
  facet_wrap(~ subject.id) +
  labs(x = "PC1", y = "PC2") +
  scale_color_identity() +
  coord_equal() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

plot_pc_over_time <- pca_data %>% 
  pivot_longer(c(PC1, PC2), names_to = "pc") %>% 
  ggplot(aes(x = day, y = value, group = pc, linetype = pc, shape = pc)) +
  geom_hline(yintercept = 0, color = "grey65") +
  geom_line() +
  geom_point() +
  facet_wrap(~ subject.id) +
  scale_color_lei(palette = "five") +
  scale_x_time_squish() +
  labs(x = "Time (days)", y = "Principal components", linetype = "Dimension", shape = "Dimension") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

#plot loadings
table_classes <- sort(table(loadings$class))
small_classes <- names(table_classes[table_classes < 3])

plot_abclass_PC <- loadings %>% 
  filter(!class %in% small_classes) %>% 
  mutate(class = gsub("-", "- ", class)) %>% 
  mutate(class = gsub("phosphatidyl", "phosphatidyl ", class)) %>% 
  arrange(desc(weights)) %>% 
  ggplot() +
  ggforce::geom_circle(aes(x0 = x0, y0 = y0, r = r), alpha = 0.3, color = NA, fill = "grey95",
                       data = data.frame(x0 = 0, y0 = 0, r = 0.5)) +
  geom_segment(aes(xend = PC1, yend = PC2, color = class, x = 0, y = 0), 
               arrow = arrow(length = unit(0.03, "npc")), show.legend = F) +
  scale_color_lei(palette = "nine", discrete = T) +
  coord_equal() +
  labs(x = "PC1", y = "PC2") +
  facet_wrap(~class, labeller = labeller(class = label_wrap_gen(15))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) 

plot_weights_pca <- loadings %>% 
  filter(!class %in% small_classes) %>% 
  pivot_longer(c(PC1, PC2), names_to = "PC") %>% 
  ggplot(aes(y = class, x = abs(value), color = class)) +
  geom_boxplot(show.legend = F) +
  geom_vline(xintercept = 0.5, color = "grey65", alpha = 0.3) +
  scale_fill_lei(palette = "nine", discrete = T) +
  scale_color_lei(palette = "nine", discrete = T) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(limits=rev) +
  labs(x = "Importance in PCA (absolute loadings)", y = "Biochemical class") +
  facet_wrap(~ PC) +
  theme_bw() +
  theme(plot.margin = unit(c(5.5, 15.5, 5.5, 5.5), "points"))+
  theme(strip.background = element_rect(fill = "white")) 

plot_weights_pca_alt <- loadings %>% 
  filter(!class %in% small_classes) %>% 
  pivot_longer(c(PC1, PC2), names_to = "PC") %>% 
  ggplot(aes(y = class, x = (value)^2*100, color = class)) +
  geom_boxplot(show.legend = F) +
  geom_vline(xintercept = 0.5^2*100, color = "grey65", alpha = 0.3) +
  scale_fill_lei(palette = "nine", discrete = T) +
  scale_color_lei(palette = "nine", discrete = T) +
  scale_x_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(limits=rev) +
  labs(x = "Importance in PCA (variance accounted for)", y = "Biochemical class") +
  facet_wrap(~ PC) +
  theme_bw() +
  theme(plot.margin = unit(c(5.5, 15.5, 5.5, 5.5), "points"))+
  theme(strip.background = element_rect(fill = "white")) 


#Gather figures
figure1 <- plot_pca_per_timepoint
figureS1 <- plot_pca_individuals
figureS2 <- plot_pc_over_time
figure2a <- plot_weights_pca
figure2b <- plot_abclass_PC


save(figure1, figureS1, figureS2, figure2a, figure2b, file = "manuscript/figures/plots_PCA.Rdata")

