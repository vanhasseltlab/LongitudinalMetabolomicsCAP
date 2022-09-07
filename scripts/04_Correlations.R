#Correlation analysis
library(tidyverse)
library(readxl)
source("functions/AllFunctions.R")
source("functions/LeidenColoring.R")

# Load data from 01_DataPreparation (data.pretreated)
load("data/01_data_clean.Rdata")

#Scale metabolite values
data.pretreated <- data.pretreated %>% 
  rename(PCT = pct, CRP = crp, CREA = crea)

metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))
names_and_class <- read.csv("data/metabolite_names_and_classes.csv", stringsAsFactors = F)
names_and_class <- names_and_class %>% distinct()

var_of_interest <- c("curb", "hospitalization.time")


#### curb and hospitalization time vs difference d0-d1, d0-d2, d0-d4 and d0-d30, metabolite values ####
data.eng <- data.pretreated %>% 
  pivot_longer(c(all_of(metrange), CRP, PCT, CREA), names_to = "metabolite") %>% 
  select(metabolite, value, day, subject.id, hospitalization.time, curb) %>% 
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
         metabolite = gsub("day\\d+_\\d+_", "",  time_metabolite)) %>% 
  select(-time_metabolite) %>% 
  pivot_wider(names_from = time, values_from = c("curb", "hospitalization.time"))

#####

#### correlations with PCT, CRP and Creatinine over time ####
# Select columns of interest
var_tdp <- c("CREA", "CRP", "PCT")
metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))

dat_met_tdp <- data.pretreated[, c("subject.id", "day", var_tdp, metrange)]
# Make correlation matrix:
cor_tdp <- cor(dat_met_tdp[, c(metrange, var_tdp)], dat_met_tdp[, var_tdp], 
               use = "pairwise.complete.obs", method = "pearson") %>% 
  as.data.frame() %>% 
  rownames_to_column("metabolite")

### Combine results #####
cor_all <- cor_tdp %>% left_join(cor_mat2)


####Visualiza curb correlations####
selected_metabolites <- c(slice_min(cor_all, curb_day30, n = 6)$metabolite)
cor_all_names <- cor_all %>% 
  left_join(names_and_class)
facet_names <- paste0(cor_all_names$name, " \ntau = ", round(cor_all_names$curb_day30, 3))
names(facet_names) <- cor_all_names$metabolite

boxplot_curb_score_selection <- data.pretreated %>%
  distinct(subject.id, .keep_all = T) %>% 
  pivot_longer(c(all_of(metrange), CRP, PCT, CREA), names_to = "metabolite") %>% 
  mutate(metabolite = factor(metabolite, levels = selected_metabolites)) %>% 
  filter(metabolite %in% selected_metabolites) %>% 
  ggplot(aes(x = curb, y = value, group = metabolite)) +
  geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names),  nrow = 2) +
  labs(y = "Scaled metabolite level", x = "CURB score") +
  theme_bw()

#####
###Heatmap only PCT and CRP
plot_correlation_heatmap <- function(correlation_matrix, 
                                     variables = names(correlation_matrix)[-1], 
                                     var_tdp = NULL, add = NULL, just_data = FALSE) {
  heatmap_dat_3 <- correlation_matrix %>% 
    pivot_longer(-metabolite, names_to = "marker", values_to = "correlation") %>% 
    filter(marker %in% variables) %>% 
    filter(metabolite %in% c(setdiff(metabolite[abs(correlation) > 0.55], variables), add)) %>% 
    left_join(names_and_class) %>% 
    mutate(name = ifelse(is.na(name), metabolite, name)) %>% 
    mutate(marker = factor(marker, levels = variables, 
                           labels = gsub("_", " - ", gsub(".", " ", variables, fixed = T), fixed = T))) %>% 
    mutate(group = ifelse(as.character(marker) %in% var_tdp, "Biomarker", gsub("[\\ \\-\\ ].*", "", marker))) %>% 
    mutate(hor_facet = ifelse(name %in% add, "1", "0")) %>% 
    distinct()
  
  m <- as.matrix(heatmap_dat_3 %>% select(name, marker, correlation) %>% 
                   pivot_wider(names_from = name, values_from = correlation) %>% 
                   select(-marker))
  clust <- hclust(dist(t(m)))
  
  if (just_data) {
    return(list(data = heatmap_dat_3, m = m, clust = clust))
  }
  
  corposterplot <- heatmap_dat_3 %>% 
    ggplot(aes(y = name, x = marker)) +
    geom_tile(aes(fill = as.numeric(correlation)))+
    scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#F54C00",
                         limits = c(-1, 1), 
                         #breaks = seq(-0.75, 0.75, 0.25),
                         breaks = seq(-1, 1, 0.25),
                         labels = seq(-1, 1, 0.25))+
    scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult  = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult  = c(0, 0)), limits= colnames(m)[clust$order]) +
    
    
    facet_grid(. ~ group, scales = "free", space='free') +
    labs(x = NULL, y = "Metabolite") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank())
  
  

  if (is.null(add)) {
    return(corposterplot)
  } else {
    intercepts <- sort(which(colnames(m)[clust$order] %in% add))
    return(corposterplot + geom_hline(yintercept = c(intercepts - 0.5, intercepts + 0.5)) )
  }
}

# Plot correlations:
#pdf("results/figures/CRP_metabolites_top10.pdf", width = 10)
print(dat_met_tdp %>% 
        pivot_longer(-c(subject.id, day, CRP), names_to = "metabolite") %>%
        filter(metabolite %in% cor_all$metabolite[abs(cor_all$CRP) > 0.6]) %>% 
        #mutate(metabolite = factor(metabolite, levels = top_and_bottom)) %>% 
        ggplot(aes(x = CRP, y = value, group = metabolite)) +
        geom_smooth(aes(x = CRP), se = FALSE, method = "lm", alpha = 0.5) +
        geom_point() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        labs(x = "C-reactive protein (CRP)", y = "Metabolite level", size = 20) +
        facet_wrap(~ metabolite, scales = "free_y") +
        theme_bw())
#dev.off()
dat_cor_plot <- plot_correlation_heatmap(cor_all, variables = grep("hospitalization.time", names(cor_all), value = T), 
                                         just_data = T, add = c("CRP", "PCT"))
intercepts <- sort(which(colnames(dat_cor_plot$m)[dat_cor_plot$clust$order] %in% c("CRP", "PCT")))

hosp_metabolite_cor <- dat_cor_plot$data %>% 
  mutate(day = gsub("hospitalization time - day", "Day ", marker, fixed = T)) %>% 
  mutate(day = factor(day, levels = paste("Day", c(1, 2, 4, 30)))) %>% 
  ggplot(aes(y = name, x = group)) +
  geom_tile(aes(fill = as.numeric(correlation)))+
  geom_hline(yintercept = c(intercepts - 0.5, intercepts + 0.5)) +
  scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#F54C00",
                       limits = c(-1, 1), 
                       #breaks = seq(-0.75, 0.75, 0.25),
                       breaks = seq(-1, 1, 0.25),
                       labels = seq(-1, 1, 0.25))+
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult  = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult  = c(0, 0)), limits= colnames(dat_cor_plot$m)[dat_cor_plot$clust$order]) +
  facet_grid(. ~ day, scales = "free", space='free') +
  labs(x = "Length of stay", y = "Metabolite") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#select higher correlations hosp time
highest_hosp_cor <- cor_all %>% 
  filter(abs(hospitalization.time_day2) > 0.551) %>% 
  arrange(desc(abs(hospitalization.time_day2))) %>% 
  select(metabolite) %>% unlist %>% as.character

data_figure4b <- data.pretreated %>% 
  select(all_of(highest_hosp_cor), subject.id, day, hospitalization.time) %>% 
  pivot_longer(all_of(highest_hosp_cor), names_to = "metabolite", values_to = "value") %>% 
  pivot_wider(names_from = day, values_from = value) %>% 
  mutate(`1` = `1` - `0`,
         `2` = `2` - `0`,
         `4` = `4` - `0`,
         `30` = `30` - `0`,
         `0` = 0) %>% 
  pivot_longer(all_of(as.character(c(0, 1, 2, 4, 30))), names_to = "day")

data_figure4b %>% 
  mutate(day = as.numeric(day)) %>% 
  #filter(day != 30) %>% 
  #filter(day %in% c(0, 2)) %>% 
  ggplot(aes(x = day, y = value, color = hospitalization.time)) +
  geom_point() +
  geom_line(aes(group = subject.id)) +
  scale_color_lei(discrete = F, palette = "gradient2") +
  scale_x_time_squish() +
  facet_grid(~ metabolite) +
  theme_bw()

rownames(names_and_class) <- names_and_class$metabolite
highest_cor_hosp_time <- data_figure4b %>% 
  left_join(names_and_class) %>%
  mutate(name = factor(name, levels = names_and_class[highest_hosp_cor, "name"] )) %>% 
  filter(day == 2) %>% 
  ggplot(aes(x = hospitalization.time, y = value)) +
  geom_smooth(method = "lm", se = F, color = "grey55", size = 0.8) +
  geom_point(aes(color = subject.id), show.legend = F) +
  labs(x = "Length of stay", y = "Change in metabolite value from day 0 to day 2") +
  scale_color_lei(discrete = T, palette = "nine") +
  facet_wrap(~ name, nrow = 2) +
  theme_bw()

#Gather figures
figure5a <- plot_correlation_heatmap(cor_all, variables = c("CRP", "PCT"), var_tdp = c("CRP", "PCT"),)
figureS2 <- plot_correlation_heatmap(cor_all, var_tdp = c("CREA", "CRP", "PCT"), add = c("CRP", "PCT", "CREA"))
figure7a <- hosp_metabolite_cor
figure7b <- highest_cor_hosp_time
figure6 <- boxplot_curb_score_selection

save(figure5a, figureS2, figure6, figure7a, figure7b, file = "manuscript/figures/plots_correlations.Rdata")

