#Correlation analysis
library(tidyverse)
source("functions/AllFunctions.R")

## Input data -----------------------------------------------------------------
# Load metabolomics data after quality control
dat.raw <- read.csv2("data/00_data_raw.csv")
metabolite_names <- read_excel("data/Metabolite_names_V2_incl_ratios.xlsx") %>% 
  select(Detected_metabolite_name_in_R, Metabolite_name_V2) %>% 
  dplyr::rename(metabolite = Detected_metabolite_name_in_R, metabolite_name = Metabolite_name_V2) %>% 
  add_row(metabolite = "LPC.20.3.", metabolite_name = "LPC (20:3)") %>% 
  mutate(metabolite_name = ifelse(metabolite_name == "R-3-Hydroxyisobutyric acid", "3-Hydroxyisobutyric acid", metabolite_name))
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
data.ratios <- AddRatiosAndSums(AddCRPAndPCT(data.reduced))
metrange <- attr(data.reduced, "metrange")




#### curb and hospitalization time vs difference d0-d30 metabolite values ####
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

#####

#### curb and hospitalization time vs difference d0-d1, d0-d2, d0-d4 and d0-d30, metabolite values ####
#TODO: add sums and ratios?
data.eng <- data.ratios %>% 
  pivot_longer(c(all_of(metrange), crp, pct, crea), names_to = "metabolite") %>% 
  select(metabolite, value, day, subject.id, hospitalization.time, curb) %>% 
  mutate(metabolite = gsub("_", ".", metabolite)) %>% 
  pivot_wider(names_from = day, values_from = value, names_prefix = "day") %>% 
  mutate(day1_0 = day1 - day0, day2_0 = day2 - day0, day4_0 = day4 - day0, day30_0 = day30 - day0) %>% 
  select(-all_of(paste0("day", c(0, 1, 2, 4, 30)))) %>% 
  pivot_wider(names_from = metabolite, values_from = c(day1_0, day2_0, day4_0, day30_0))

metrange_new <- setdiff(names(data.eng), c("subject.id", "day", var_of_interest))

cor_curb <- cor(data.eng[, metrange_new], data.eng[, "curb"], use = "pairwise.complete.obs", method = "kendall")
cor_ht <- cor(data.eng[, metrange_new], data.eng[, "hospitalization.time"], use = "pairwise.complete.obs", method = "pearson")
cor_mat2 <- as.data.frame(cbind(cor_curb, cor_ht)) %>% rownames_to_column("metabolite") %>% 
  dplyr::rename(time_metabolite = metabolite) %>% 
  mutate(time = str_match(time_metabolite, "day\\d+")[, 1],
         metabolite = str_match(time_metabolite, "([^\\_]+$)")[, 1]) %>% 
  select(-time_metabolite) %>% 
  pivot_wider(names_from = time, values_from = c("curb", "hospitalization.time"))

#####

#### correlations with PCT, CRP and Creatinine over time ####
# Select columns of interest
var_tdp <- c("crea", "crp", "pct")
metrange <- c(attr(data.ratios, "metrange"), attr(data.ratios, "ratiorange"))

dat_met_tdp <- data.ratios[, c("subject.id", "day", var_tdp, metrange)]
# Make correlation matrix:
cor_tdp <- cor(dat_met_tdp[, c(metrange, var_tdp)], dat_met_tdp[, var_tdp], 
               use = "pairwise.complete.obs", method = "pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column("metabolite") %>% 
  mutate(metabolite = gsub("_", ".", metabolite))

#####

### Combine results #####
cor_all <- cor_tdp %>% left_join(cor_mat2)



#### Visualize heatmap with correlations ####
heatmap_dat_3 <- cor_all %>% 
  pivot_longer(-metabolite, names_to = "marker", values_to = "correlation") %>% 
  #filter(marker %in% c("pct", "crp")) %>% 
  filter(metabolite %in% metabolite[abs(correlation) > 0.55]) %>% 
  left_join(metabolite_names %>% mutate(metabolite = gsub("_", ".", metabolite))) %>% 
  mutate(metabolite_name = ifelse(is.na(metabolite_name), metabolite, metabolite_name)) %>% 
  mutate(marker = factor(marker, levels = names(cor_all)[-1], 
                         labels = gsub("_", " - ", gsub(".", " ", R.utils::capitalize(names(cor_all)[-1]), fixed = T), fixed = T))) %>% 
  mutate(group = ifelse(as.character(marker) %in% R.utils::capitalize(var_tdp), "Biomarker", gsub("[\\ \\-\\ ].*", "", marker)))
  
cor_selected <- cor_all %>% 
  filter(metabolite %in% heatmap_dat_3$metabolite) %>% 
  left_join(heatmap_dat_3 %>% select(metabolite, metabolite_name) %>% distinct()) 

m <- as.matrix(cor_selected %>%
                 column_to_rownames("metabolite_name") %>% 
                 select(-c(metabolite))) #%>% 
                 #select(-c(hospitalization.time_day1, hospitalization.time_day2, hospitalization.time_day4, hospitalization.time_day30))
clust <- hclust(dist(m))

corposterplot <- heatmap_dat_3 %>% 
  #filter(marker != "Crea") %>% 
  #filter(!metabolite %in% var_tdp) %>% 
  ggplot(aes(y = metabolite_name, x = marker)) +
  geom_tile(aes(fill = as.numeric(correlation)))+
  scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#FF9933",
                       limits = c(-1, 1), 
                       #breaks = seq(-0.75, 0.75, 0.25),
                       breaks = seq(-1, 1, 0.25),
                       labels = seq(-1, 1, 0.25))+
  scale_y_discrete(limits = rownames(m)[clust$order], drop = T)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(~ group, scales = "free_x", space='free') +
  labs(x = NULL, y = "Metabolite") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.x = element_blank())


pdf(file = "results/figures/correlations_markers_metabolites.pdf", height = 10)
print(corposterplot)
dev.off()


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





