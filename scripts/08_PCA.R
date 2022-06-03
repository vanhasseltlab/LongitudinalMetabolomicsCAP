##PCA
source("functions/AllFunctions.R")
source("functions/LeidenColoring.R")


# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
#Remove metabolites with to many NAs
data.clean <- DataCleaning(data.reduced, attr(data.reduced, "metrange"))

#Add longitudinal markers: CRP and PCT
data.full <- AddCRPAndPCT(data.reduced)

#Add ratios and sums to data.full
data.full.clean <- AddRatiosAndSums(AddCRPAndPCT(data.clean))

#Log transform and scale metabolite data
data.full.scaled <- DataPretreatment(data.full.clean, metrange = attr(data.clean, "metrange"), scaling = "auto")


metrange <- attr(data.full.scaled, "metrange")
metrange <- metrange[colSums(is.na(data.full.scaled[, metrange])) < 2]

pca_metabolites <- data.full.scaled[!is.na(rowSums(data.full.scaled[, metrange])), metrange]
set.seed(123)
pca_all_times <- prcomp(pca_metabolites)

pca_data <- data.frame(data.full.scaled[rownames(pca_metabolites), c("day", "subject.id", "hospitalization.time", "curb", "age", "sex")], 
                       pca_all_times$x[, 1:4])

plot_pca_per_timepoint <- pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = hospitalization.time)) +
  geom_point() +
  scale_color_viridis_c(direction = -1) +
  facet_grid(~ day) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom")

plot_pca_per_timepoint_curb <- pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(curb))) +
  geom_point() +
  scale_color_viridis_d(direction = -1) +
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
  #geom_point(data = pca_data, aes(x = PC1, y = PC2), show.legend = F) +

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
  #geom_point(data = pca_data, aes(x = PC1, y = PC2), show.legend = F) +
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = day), show.legend = F) +
  geom_segment(aes(x = PC1_0, y = PC2_0, xend = PC1_1, yend = PC2_1)) +
  geom_segment(aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2)) +
  geom_segment(aes(x = PC1_2, y = PC2_2, xend = PC1_4, yend = PC2_4)) +
  geom_segment(aes(x = PC1_4, y = PC2_4, xend = PC1_30, yend = PC2_30)) +
  theme_bw()

pca_plot_hosp_time <- pca_data %>% 
  pivot_wider(names_from = day, values_from = c(PC1, PC2, PC3, PC4)) %>% 
  ggplot(aes(group = subject.id, color = hospitalization.time)) +
  #geom_point(data = pca_data, aes(x = PC1, y = PC2), show.legend = F) +
  
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


pdf(file = "results/figures/PCA_time_profiles.pdf", height = 10, width = 10)
print(pca_plot)
dev.off()

pca_plot
ggsave("results/figures/PCA_time_profiles.png", width = 10, height = 10, dpi = 400)

pdf(file = "results/figures/PCA_exploration.pdf", height = 5, width = 7)
print(pca_plot_hosp_time)
print(pca_plot_curb)
print(pca_plot_age)
print(pca_plot_sex)
print(pca_plot_facetted)
print(plot_PC_time_curve("PC1"))
print(plot_PC_time_curve("PC2"))
print(plot_PC_time_curve("PC3"))
print(plot_PC_time_curve("PC4"))
dev.off()



#Gather figures
figure6a <- plot_pca_per_timepoint
figure6b <- plot_PC_time_curve("PC1", color = "hospitalization.time") +
  scale_color_viridis_c(direction = -1)
figure6c <- plot_PC_time_curve("PC2", color = "hospitalization.time") +
  scale_color_viridis_c(direction = -1)
figureS6 <- pca_plot

save(figure6a, figure6b, figure6c, figureS6, file = "manuscript/figures/plots_PCA.Rdata")
