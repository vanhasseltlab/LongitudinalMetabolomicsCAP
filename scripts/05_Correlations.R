#correlation analysis
library(tidyverse)

#data.raw <- read.csv2("data/00_data_reduced.csv")
load("data/column_names.Rdata")

#data.reduced from 10_generalscript.R
data.reduced

#PSI vs difference d0-d30 metaboliete values:
var_of_interest <- c("psi.score", "curb", "duration.of.symptoms.before.admission")

dat_psi_met <- data.reduced[, c("subject.id", "day", var_of_interest, metrange)] %>% 
  filter(day %in% c(0, 30)) %>% 
  pivot_longer(-c(subject.id, day, all_of(var_of_interest)), names_to = "metabolite") %>% 
  pivot_wider(id_cols = c(subject.id, metabolite, all_of(var_of_interest)), names_from = day) %>% 
  mutate(day30_0_diff = `30` - `0`) %>% 
  pivot_wider(names_from = metabolite, values_from = day30_0_diff, id_cols = c(subject.id,  all_of(var_of_interest))) %>% 
  mutate(psi.score = factor(psi.score, levels = levels(psi.score)[order(as.numeric(gsub("-\\d*", "", levels(psi.score))))]))

cor_mat <- cor(dat_psi_met[, metrange], as.numeric(dat_psi_met$curb), use = "pairwise.complete.obs", method = "kendall")

nr_top_and_bottom <- 10
top_and_bottom <- c(as.data.frame(cor_mat) %>% slice_max(V1, n = nr_top_and_bottom) %>% rownames,
                    as.data.frame(cor_mat) %>% slice_min(V1, n = nr_top_and_bottom) %>% rownames)

lowest_cor <- c(as.data.frame(cor_mat) %>% filter(V1 < -0.5) %>% 
                  rownames_to_column('metabolite') %>% 
                  arrange(V1) %>% select(metabolite) %>% unlist())

facet_names <- paste0(rownames(cor_mat), " (tau = ", round(cor_mat[, 1], 3), ")")
names(facet_names) <- rownames(cor_mat)

pdf("results/figures/Curb_metabolite0-30_top10.pdf", width = 10)
print(dat_psi_met %>% 
  pivot_longer(-c(subject.id, all_of(var_of_interest)), names_to = "metabolite") %>%
  filter(metabolite %in% lowest_cor) %>% 
  mutate(metabolite = factor(metabolite, levels = lowest_cor)) %>% 
  ggplot(aes(x = curb, y = value, group = metabolite)) +
  geom_smooth(aes(x = curb), se = FALSE, method = "lm", alpha = 0.5) +
  geom_point() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "PSI score", y = NULL) +
  facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
  theme_bw())
dev.off()

pdf("results/figures/Curb_metabolite0-30_top10.pdf", width = 10)
print(dat_psi_met %>% 
        pivot_longer(-c(subject.id, all_of(var_of_interest)), names_to = "metabolite") %>%
        filter(metabolite %in% lowest_cor) %>% 
        filter(!is.na(value)) %>% 
        mutate(metabolite = factor(metabolite, levels = lowest_cor)) %>% 
        ggplot(aes(x = curb, y = value, group = metabolite)) +
        geom_boxplot(aes(x = as.factor(curb), group =  as.factor(curb))) +
        labs(x = "Curb score", y = NULL) +
        facet_wrap(~ metabolite, scales = "free_y", labeller = labeller(metabolite = facet_names)) +
        theme_bw())
dev.off()

hist(cor_mat[, 1], freq = FALSE, breaks = 30)
x <- seq(-1, 1, by = 0.05)
lines(x, dnorm(x, 0, sd(cor_mat[, 1])), col = "red")

