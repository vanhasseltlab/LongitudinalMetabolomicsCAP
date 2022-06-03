#Correlation analysis
library(tidyverse)
library(readxl)
source("functions/AllFunctions.R")

## Input data -----------------------------------------------------------------
dat.raw <- read.csv2("data/00_data_raw.csv")

#Data preparation
# Remove redundant variables from data frame
data.reduced <- ReduceData(dat.raw)
#Remove metabolites with to many NAs
data.clean <- DataCleaning(data.reduced, attr(data.reduced, "metrange"))

#Add longitudinal markers: CRP and PCT and ratios and sums
data.pc <- AddCRPAndPCT(data.clean)
data.rs <- AddRatiosAndSums(data.pc)

#Scale metabolite values
data.pretreated <- DataPretreatment(data.rs, metrange = c(attr(data.rs, "metrange"), attr(data.rs, "ratiorange"))) %>% 
  rename(PCT = pct, CRP = crp, CREA = crea)

metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"))
names_and_class <- read.csv("data/metabolite_names_and_classes.csv", stringsAsFactors = F)

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
selected_metabolites <- c(slice_max(cor_all, curb_day30, n = 4)$metabolite, 
                          slice_min(cor_all, curb_day30, n = 4)$metabolite)
cor_all_names <- cor_all %>% 
  left_join(names_and_class)
facet_names <- paste0(cor_all_names$name, " \ntau = ", round(cor_all$curb_day30, 3))
names(facet_names) <- cor_all_names$metabolite

boxplot_curb_score_selection <- data.pretreated %>%
  distinct(subject.id, .keep_all = T) %>% 
  pivot_longer(c(all_of(metrange), CRP, PCT, CREA), names_to = "metabolite") %>% 
  mutate(metabolite = factor(metabolite, levels = selected_metabolites)) %>% 
  filter(metabolite %in% selected_metabolites) %>% 
  ggplot(aes(x = curb, y = value, group = metabolite)) +
  geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names),  nrow = 2) +
  theme_bw()

#####
###Heatmap only PCT and CRP
plot_correlation_heatmap <- function(correlation_matrix, variables = names(correlation_matrix)[-1], var_tdp = NULL, add = NULL) {
  heatmap_dat_3 <- correlation_matrix %>% 
    pivot_longer(-metabolite, names_to = "marker", values_to = "correlation") %>% 
    filter(marker %in% variables) %>% 
    filter(metabolite %in% c(setdiff(metabolite[abs(correlation) > 0.55], variables), add)) %>% 
    left_join(names_and_class) %>% 
    mutate(name = ifelse(is.na(name), metabolite, name)) %>% 
    mutate(marker = factor(marker, levels = variables, 
                           labels = gsub("_", " - ", gsub(".", " ", variables, fixed = T), fixed = T))) %>% 
    mutate(group = ifelse(as.character(marker) %in% var_tdp, "Biomarker", gsub("[\\ \\-\\ ].*", "", marker))) %>% 
    mutate(hor_facet = ifelse(name %in% add, "1", "0"))
  
  corposterplot <- heatmap_dat_3 %>% 
    ggplot(aes(y = name, x = marker)) +
    geom_tile(aes(fill = as.numeric(correlation)))+
    scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#FF9933",
                         limits = c(-1, 1), 
                         #breaks = seq(-0.75, 0.75, 0.25),
                         breaks = seq(-1, 1, 0.25),
                         labels = seq(-1, 1, 0.25))+
    scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult  = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult  = c(0, 0))) +
    facet_grid(hor_facet ~ group, scales = "free", space='free') +
    labs(x = NULL, y = "Metabolite") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank())
  
  return(corposterplot)
}

plot_correlation_heatmap(cor_all, var_tdp = c("CREA", "CRP", "PCT"))
plot_correlation_heatmap(cor_all, variables = c("CRP", "PCT"), var_tdp = c("CRP", "PCT"))
plot_correlation_heatmap(cor_all, variables = grep("hospitalization.time", names(cor_all), value = T), add = c("CRP", "PCT"))


pdf(file = "results/figures/correlations_markers_metabolites.pdf", height = 10)
plot_correlation_heatmap(cor_all, variables = grep("hospitalization.time", names(cor_all), value = T), add = c("CRP", "PCT"))
dev.off()


# Add correlations to general dataframe
#cor_mat_2 <-data.frame(metabolite = rownames(cor_mat), marker = "CRP", correlation = as.data.frame(cor_mat)$V1)

# Make subselection of data for plotting:
# nr_top_and_bottom <- 3
# top_and_bottom <- c(as.data.frame(cor_mat) %>% slice_max(V1, n = nr_top_and_bottom) %>% rownames,
#                     as.data.frame(cor_mat) %>% slice_min(V1, n = nr_top_and_bottom) %>% rownames)
# Plot correlations:
facet_names <- paste0(rownames(cor_mat), " (cor = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)
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



#Gather figures
figure3 <- plot_correlation_heatmap(cor_all, variables = c("CRP", "PCT"), var_tdp = c("CRP", "PCT"))
figureS4 <- plot_correlation_heatmap(cor_all, var_tdp = c("CREA", "CRP", "PCT"), add = c("CRP", "PCT", "CREA"))
figureS5 <- plot_correlation_heatmap(cor_all, variables = grep("hospitalization.time", names(cor_all), value = T), add = c("CRP", "PCT"))
figure5 <- boxplot_curb_score_selection

save(figure3, figureS4, figureS5, figure5, file = "manuscript/figures/plots_correlations.Rdata")

