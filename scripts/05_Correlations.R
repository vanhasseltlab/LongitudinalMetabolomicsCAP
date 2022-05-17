#Correlation analysis
library(tidyverse)
source("functions/AllFunctions.R")

## Input data -----------------------------------------------------------------
# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
metrange <- attr(data.reduced, "metrange")

#curb and hospitalization time vs difference d0-d30 metabolite values:
var_of_interest <- c("curb", "hospitalization.time")

dat_curb_met <- data.reduced[, c("subject.id", "day", var_of_interest, metrange)] %>% 
  filter(day %in% c(0, 30)) %>% 
  pivot_longer(-c(subject.id, day, all_of(var_of_interest)), names_to = "metabolite") %>% 
  pivot_wider(id_cols = c(subject.id, metabolite, all_of(var_of_interest)), names_from = day) %>% 
  mutate(day30_0_diff = `30` - `0`) %>% 
  pivot_wider(names_from = metabolite, values_from = day30_0_diff, id_cols = c(subject.id,  all_of(var_of_interest)))

cor_curb <- cor(dat_curb_met[, metrange], dat_curb_met[, "curb"], use = "pairwise.complete.obs", method = "kendall")
cor_ht <- cor(dat_curb_met[, metrange], dat_curb_met[, "hospitalization.time"], use = "pairwise.complete.obs", method = "pearson")
cor_mat <- as.data.frame(cbind(cor_curb, cor_ht)) %>% rownames_to_column("metabolite")

lowest_cor <- cor_mat %>% filter(abs(curb) > 0.55 | abs(hospitalization.time) > 0.55) %>%  
                  arrange(curb) %>% select(metabolite) %>% unlist()

facet_names <- paste0(cor_mat$metabolite, " \nCURB tau = ", round(cor_mat$curb, 3), ",\nhosp rho = ", round(cor_mat$hospitalization.time, 3))
names(facet_names) <- cor_mat$metabolite

plot_data <- dat_curb_met %>% 
  pivot_longer(-c(subject.id, all_of(var_of_interest)), names_to = "metabolite") %>%
  filter(metabolite %in% lowest_cor) %>% 
  filter(!is.na(value)) %>% 
  mutate(metabolite = factor(metabolite, levels = lowest_cor)) 

pdf("results/figures/Curb_metabolite0-30_top10.pdf", width = 10)
print(plot_data %>% 
  ggplot(aes(x = hospitalization.time, y = value, group = metabolite)) +
  geom_point(aes(color = as.factor(curb))) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
  theme_bw())
print(plot_data %>% 
  ggplot(aes(x = curb, y = value, group = metabolite)) +
  geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
  theme_bw())
dev.off()


#correlations with PCT, CRP and Creatinine over time
# Select columns of interest
data.ratios <- AddRatiosAndSums(data.full)
var_tdp <- c("crp", "crea", "pct")
metrange <- c(attr(data.ratios, "metrange"), attr(data.ratios, "ratiorange"))

dat_met_tdp <- data.ratios[, c("subject.id", "day", var_tdp, metrange)]
# Make correlation matrix:
cor_tdp <- cor(dat_met_tdp[, metrange], dat_met_tdp[, var_tdp], use = "pairwise.complete.obs", method = "pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column("metabolite")



###Recreate poster heatmap
metabolite_names <- read_excel("data/Metabolite_names_V2_incl_ratios.xlsx") %>% 
  select(Detected_metabolite_name_in_R, Metabolite_name_V2) %>% 
  dplyr::rename(metabolite = Detected_metabolite_name_in_R, metabolite_name = Metabolite_name_V2) %>% 
  add_row(metabolite = "LPC.20.3.", metabolite_name = "LPC (20:3)") %>% 
  mutate(metabolite_name = ifelse(metabolite_name == "R-3-Hydroxyisobutyric acid", "3-Hydroxyisobutyric acid", metabolite_name))

heatmap_dat_3_test <- cor_tdp %>% 
  pivot_longer(-metabolite, names_to = "marker", values_to = "correlation") %>% 
  #filter(marker %in% c("pct", "crp")) %>% 
  filter(metabolite %in% metabolite[abs(correlation) > 0.55]) %>% 
  left_join(metabolite_names) %>% 
  mutate(metabolite_name = ifelse(is.na(metabolite_name), metabolite, metabolite_name))
  
cor_selected <- cor_tdp %>% 
  filter(metabolite %in% heatmap_dat_3_test$metabolite) %>% 
  left_join(heatmap_dat_3_test %>% select(metabolite, metabolite_name) %>% distinct())
m <- as.matrix(cor_selected[, c("pct", "crp", "crea")])
clust <- hclust(dist(m))

corposterplot <- ggplot(data = heatmap_dat_3_test, aes(y = metabolite_name, x = marker)) +
  geom_tile(aes(fill = as.numeric(correlation)))+
  scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#FF9933",
                       limits = c(-0.75, 0.75), 
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75),
                       labels = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75))+
  scale_y_discrete(limits = cor_selected[clust$order, "metabolite_name"])+
  theme_bw()+
  theme(text=element_text(size=20))

# Add correlations to general dataframe
#cor_mat_2 <-data.frame(metabolite = rownames(cor_mat), marker = "CRP", correlation = as.data.frame(cor_mat)$V1)
corr_all <- cor_tdp %>% 
  left_join(cor_mat)
# Make subselection of data for plotting:
# nr_top_and_bottom <- 3
# top_and_bottom <- c(as.data.frame(cor_mat) %>% slice_max(V1, n = nr_top_and_bottom) %>% rownames,
#                     as.data.frame(cor_mat) %>% slice_min(V1, n = nr_top_and_bottom) %>% rownames)
# Plot correlations:
facet_names <- paste0(rownames(cor_mat), " (cor = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)
#pdf("results/figures/CRP_metabolites_top10.pdf", width = 10)
print(dat_met_tdp %>% 
        pivot_longer(-c(subject.id, day, crp), names_to = "metabolite") %>%
        filter(metabolite %in% corr_all$metabolite[abs(corr_all$crp) > 0.6]) %>% 
        #mutate(metabolite = factor(metabolite, levels = top_and_bottom)) %>% 
        ggplot(aes(x = crp, y = value, group = metabolite)) +
        geom_smooth(aes(x = crp), se = FALSE, method = "lm", alpha = 0.5) +
        geom_point() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        labs(x = "C-reactive protein (CRP)", y = "Metabolite level", size = 20) +
        facet_wrap(~ metabolite, scales = "free_y") +
        theme_bw())
#dev.off()





