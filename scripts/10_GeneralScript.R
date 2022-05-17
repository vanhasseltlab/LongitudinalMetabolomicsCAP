## With this script, all analysis for the Longitunal Metabolomics Project, can be reproduced
## Code developped by Ilona den Hartog and Laura Zwep

## Setup R --------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Load libraries
library(plyr)
library(foreign)
library(ggplot2)
library(openxlsx)
library(tidyverse)
library(readxl)
# Source functions
source("functions/AllFunctions.R")

## Input data -----------------------------------------------------------------
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


## Load metabolite class data


## Define variables ----

# Define metabolite and metadata range in data after quality control
metrange <- attr(data.full, "metrange")
nonmetrange <- setdiff(names(data.full), metrange)

## Correlations ----

## Change in metabolite level (d0-d30) related to disease severity:
# Change in metabolite levels (d0-d30) vs CURB score (d0)
# Calulate change metabolite levels (d30-d0)
dat_met_d30_d0 <- data.full[, c("subject.id", "day", metrange)] %>% 
  filter(day %in% c(0, 30)) %>% 
  pivot_longer(-c(subject.id, day), names_to = "metabolite") %>% 
  pivot_wider(id_cols = c(subject.id, metabolite), names_from = day) %>% 
  mutate(met_d30_d0 = `30` - `0`) %>% 
  pivot_wider(names_from = metabolite, values_from = met_d30_d0, id_cols = c(subject.id))

# Select 2nd variable for correlation calculation
dat_curb <- select(data.full, "subject.id", "day", "curb") %>% 
  filter(day == 0) %>% 
  select(-day)

# Combine variables of interst to dataframe for correlation calculation
cor_dat <- left_join(dat_curb, dat_met_d30_d0, by = c("subject.id"))

# Make correlation matrix:
cor_mat <- cor(cor_dat[, metrange], as.numeric(cor_dat$curb), use = "pairwise.complete.obs", method = "kendall") # method kendall for ordinal variable CURB
# add results to general correlation dataframe
cor_mat_2 <-data.frame(metabolite = rownames(cor_mat), marker = "CURB", correlation = as.data.frame(cor_mat)$V1)
corr_all <- cor_mat_2

#metabolite names
metabolite_names <- read_excel("data/Metabolite_names_V2.xlsx") %>% 
  as.data.frame(metabolite_names) %>% 
  select(Detected_metabolite_name_in_R, Metabolite_name_V2)
colnames(metabolite_names) <- c("metabolite", "metabolite_name")
metabolite_names$metabolite_name[metabolite_names$metabolite_name == "R-3-Hydroxyisobutyric acid"] <- "3-Hydroxyisobutyric acid"


# Make subselection of data for plotting:
lowest_cor <- c(as.data.frame(cor_mat) %>% filter(V1 < -0.55) %>% 
                  rownames_to_column('metabolite') %>% 
                  arrange(V1) %>% select(metabolite) %>% unlist())

# Plot correlations:
facet_names <- paste0(rownames(cor_mat), " (tau = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)

plotdat <- cor_dat %>% 
  pivot_longer(-c(subject.id, curb), names_to = "metabolite") %>%
  filter(metabolite %in% lowest_cor) %>% 
  filter(!is.na(value)) %>%  
  mutate(metabolite = factor(metabolite, levels = lowest_cor))

posterplot <- ggplot(data = plotdat, aes(x = curb, y = value, group = metabolite)) +
  geom_abline(aes(intercept = 0, slope = 0),  linetype="dotted")+
  geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
  labs(x = "CURB score", y = "Metabolite change (day 30 - day 0)") +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names), ncol = 3) +
  theme_bw() +
  theme(text=element_text(size=28))

ggsave(posterplot, file="results/figures/CURB_met_d30-d0_t0.55_boxplot.png", width=12, height=7, dpi=300)

pdf("results/figures/CURB_met_d0-d30_t0.55_boxplot.pdf", width = 10)
print(cor_dat %>% 
        pivot_longer(-c(subject.id, curb), names_to = "metabolite") %>%
        filter(metabolite %in% lowest_cor) %>% 
        filter(!is.na(value)) %>% 
        mutate(metabolite = factor(metabolite, levels = lowest_cor)) %>% 
        ggplot(aes(x = curb, y = value, group = metabolite)) +
        geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
        labs(x = "CURB score", y = "Metabolite change (day 30 - day 0)") +
        facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names), ncol = 5) +
        theme_bw())
dev.off()


## Change in metabolite levels (d0, 1, 2, 4, 30) vs CRP change (d0, 1, 2, 4, 30)
# Select columns of interest
metrange <- attr(data.full, "metrange")
dat_met_crp <- data.full[, c("subject.id", "day", "crp", metrange)]
# Make correlation matrix:
cor_mat <- cor(dat_met_crp[, metrange], as.numeric(dat_met_crp$crp), use = "pairwise.complete.obs", method = "pearson")
# Add correlations to general dataframe
cor_mat_2 <-data.frame(metabolite = rownames(cor_mat), marker = "CRP", correlation = as.data.frame(cor_mat)$V1)
corr_all <- corr_all %>% 
  bind_rows(cor_mat_2)
# Make subselection of data for plotting:
nr_top_and_bottom <- 3
top_and_bottom <- c(as.data.frame(cor_mat) %>% slice_max(V1, n = nr_top_and_bottom) %>% rownames,
                    as.data.frame(cor_mat) %>% slice_min(V1, n = nr_top_and_bottom) %>% rownames)
# Plot correlations:
facet_names <- paste0(rownames(cor_mat), " (cor = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)
pdf("results/figures/CRP_metabolites_top10.pdf", width = 10)
print(dat_met_crp %>% 
        pivot_longer(-c(subject.id, day, crp), names_to = "metabolite") %>%
        filter(metabolite %in% top_and_bottom) %>% 
        mutate(metabolite = factor(metabolite, levels = top_and_bottom)) %>% 
        ggplot(aes(x = crp, y = value, group = metabolite)) +
        geom_smooth(aes(x = crp), se = FALSE, method = "lm", alpha = 0.5) +
        geom_point() +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        labs(x = "C-reactive protein (CRP)", y = "Metabolite level", size = 20) +
        facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
        theme_bw())
dev.off()


## Change in metabolite levels (d0, 1, 2, 4, 30) vs PCT change (d0, 1, 2, 4, 30)
# Select columns of interest
dat_met_pct <- data.full[, c("subject.id", "day", "pct", metrange)]
# Make correlation matrix:
cor_mat <- cor(dat_met_pct[, metrange], as.numeric(dat_met_pct$pct), use = "pairwise.complete.obs", method = "pearson")
# Add correlations to general dataframe
cor_mat_2 <-data.frame(metabolite = rownames(cor_mat), marker = "PCT", correlation = as.data.frame(cor_mat)$V1)
corr_all <- corr_all %>% 
  bind_rows(cor_mat_2)
# Make subselection of data for plotting:
nr_top_and_bottom <- 10
top_and_bottom <- c(as.data.frame(cor_mat) %>% slice_max(V1, n = nr_top_and_bottom) %>% rownames,
                    as.data.frame(cor_mat) %>% slice_min(V1, n = nr_top_and_bottom) %>% rownames)
# Plot correlations:
facet_names <- paste0(rownames(cor_mat), " (cor = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)
pdf("results/figures/PCT_metabolites_top10.pdf", width = 10)
print(dat_met_pct %>% 
        pivot_longer(-c(subject.id, day, pct), names_to = "metabolite") %>%
        filter(metabolite %in% top_and_bottom) %>% 
        mutate(metabolite = factor(metabolite, levels = top_and_bottom)) %>% 
        ggplot(aes(x = pct, y = value, group = metabolite)) +
        geom_smooth(aes(x = pct), se = FALSE, method = "lm", alpha = 0.5) +
        geom_point() +
        # scale_x_discrete(guide = guide_axis(angle = 90)) +
        labs(x = "Procalcitonin", y = "Metabolite level") +
        facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
        theme_bw())
dev.off()


## PCT (d0, 1, 2, 4, 30) vs CRP (d0, 1, 2, 4, 30)
# Select columns of interest
dat_crp_pct <- data.full[, c("subject.id", "day", "pct", "crp")]
# Make correlation matrix:
cor_mat <- cor(as.numeric(dat_crp_pct$crp), as.numeric(dat_crp_pct$pct), use = "pairwise.complete.obs", method = "pearson")
# Correlation is 0.6 -> CRP and PCT are correlated


## CRP (d0) vs CURB (d0)
# Select columns of interest
dat_crp_curb <- select(data.full, subject.id, day, crp, curb) %>% 
  filter(day == 0)
# Make correlation matrix:
cor_mat <- cor(as.numeric(dat_crp_curb$crp), as.numeric(dat_crp_curb$curb), use = "pairwise.complete.obs", method = "kendall")
# Correlation is 0.23


## PCT (d0) vs CURB (d0)
# Select columns of interest
dat_pct_curb <- select(data.full, subject.id, day, pct, curb) %>% 
  filter(day == 0)
# Make correlation matrix:
cor_mat <- cor(as.numeric(dat_pct_curb$pct), as.numeric(dat_pct_curb$curb), use = "pairwise.complete.obs", method = "kendall")
# Correlation is 0.4

## Correlations heatmap
heatmap_dat <- corr_all %>% 
  filter(marker == "CRP" | marker == "PCT") %>% 
  filter(abs(correlation) > 0.55)
heatmap_dat_2 <- corr_all %>% 
  filter(metabolite %in% heatmap_dat$metabolite) %>% 
  filter(marker == "CRP" | marker == "PCT")


metabolite_names <- read_excel("data/Metabolite_names_V2.xlsx") %>% 
  as.data.frame(metabolite_names) %>% 
  filter(Detected_metabolite_name_in_R %in% heatmap_dat_2$metabolite) %>% 
  select(Detected_metabolite_name_in_R, Metabolite_name_V2)
colnames(metabolite_names) <- c("metabolite", "metabolite_name")
metabolite_names <- add_row(metabolite_names, metabolite = "LPC.20.3.", metabolite_name = "LPC (20:3)")
metabolite_names$metabolite_name[metabolite_names$metabolite_name == "R-3-Hydroxyisobutyric acid"] <- "3-Hydroxyisobutyric acid"

heatmap_dat_3 <- pivot_wider(heatmap_dat_2, names_from  = marker, values_from = correlation) %>% 
  left_join(metabolite_names, by = "metabolite") %>% 
  select(-metabolite) %>% 
  filter(metabolite_name != "S-3-Hydroxyisobutyric acid") %>%
  pivot_longer(names_to = "marker", values_to = "correlation", -metabolite_name) %>%
  dplyr::rename("metabolite" = "metabolite_name")

m <- as.matrix(pivot_wider(heatmap_dat_3, names_from = metabolite, values_from = correlation)[, -1])
clust <- hclust(dist(t(m)))

corposterplot <- ggplot(data = heatmap_dat_3, aes(y = metabolite, x = marker)) +
  geom_tile(aes(fill = as.numeric(correlation)))+
  scale_fill_gradient2(name = "Correlation", low = "#001158", high = "#FF9933",
                       limits = c(-0.75, 0.75), 
                       breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75),
                       labels = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75))+
  scale_y_discrete(limits = colnames(m)[clust$order])+
  scale_x_discrete(limits = c("CRP", "PCT"))+
  theme_bw()+
  theme(text=element_text(size=20))
  
ggsave(corposterplot, file="results/figures/heatmap_PCT_CRP_met_over_time.png", width=6, height=8, dpi=300)

corr_all_wide <- pivot_wider(corr_all, names_from = marker, values_from = correlation)
#Combine with spline data
spline_data_2 <- spline_data %>% 
  rownames_to_column("metabolite") %>% 
  select(metabolite, cluster)
corr_all_wide <- left_join(corr_all_wide, spline_data_2, by = "metabolite")

metabolite_names <- read_excel("data/Metabolite_names_V2.xlsx") %>% 
  as.data.frame(metabolite_names) %>% 
  filter(Detected_metabolite_name_in_R %in% corr_all_wide$metabolite) %>% 
  select(Detected_metabolite_name_in_R, Metabolite_name_V2)
colnames(metabolite_names) <- c("metabolite", "metabolite_name")

corr_all_wide <- corr_all_wide %>% 
  left_join(metabolite_names, by = "metabolite") %>% 
  select(-metabolite)
write.csv(corr_all_wide, "results/corr_all_wide.csv")
  
#Biological interpretation or CURB, CRP/PCT and splines
decreasing_trend <- corr_all_wide %>% 
  filter(cluster == 5)

