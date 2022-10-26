# Get splines for metabolomics and CRP
library(tidyverse)
library(foreign)

# Source functions
source("functions/LeidenColoring.R")

names_and_class <- read.csv("data/metabolite_names_and_classes.csv", stringsAsFactors = F)

# Load data from 01_DataPreparation (data.pretreated)
load("data/01_data_clean.Rdata")

# Load CRP data
infl.markers.longitudinal <- read.spss("data/raw/aanvullende data ilona 1762021.sav",
                                       use.value.labels = TRUE, to.data.frame = TRUE)

### Calculate splines for each metabolite, over all patients
data.pretreated$PCT <- with(data.pretreated, (pct - mean(pct, na.rm = T))/sd(pct, na.rm = T))

metrange <- c(attr(data.pretreated, "metrange"), attr(data.pretreated, "ratiorange"), "PCT")
spline_data <- as.data.frame(matrix(ncol = 5, nrow = length(metrange)))
rownames(spline_data) <- metrange
for (met in metrange) {
  dat <- data.pretreated[, c(met, "day")] %>% 
    filter(!is.na(!!as.name(met)))
  spline_met <- smooth.spline(x = dat$day, y = dat[, met])
  spline_data[met, ] <- spline_met$y
}
colnames(spline_data)[1:5] <- sort(unique(dat$day))
long_spline <- spline_data %>% 
  rownames_to_column("metabolite") %>% 
  pivot_longer(-c(metabolite), names_to = "day", values_to = "spline") %>% 
  mutate(day = as.numeric(day))

### Calculate spline for CRP
dat_crp <- select(infl.markers.longitudinal, Studienr, CRP_0, CRP_dg1, CRP_dg2, CRP_3, CRP_dg4, 
                  CRP_dg5, CRP_dg6, CRP_7, CRP_dg10, CRP_30)
dat_crp_long <- dat_crp %>% 
  pivot_longer(-c(Studienr), names_to = "day", values_to = "crp_value") %>% 
  na.omit() %>% #Remove NAs
  mutate(crp_value = log2(crp_value + 1)) %>% #Log transform
  mutate(crp_value = (crp_value - mean(crp_value)) / sd(crp_value)) %>% #Scale
  mutate(day = as.numeric(gsub("CRP_|CRP_dg", "", day)))

spline_crp <- smooth.spline(x = as.numeric(as.character(dat_crp_long$day)), y = dat_crp_long$crp_value)
long_crp <- data.frame(day = spline_crp$x, spline = spline_crp$y)

### visualize interesting metabolites
metabolites_of_interest <- c("PC.34.1.", "LPC.14.0.", "LPC.16.0.", "LPC.16.1.", "LPC.18.3.", "LPC.20.4.", "PCT")

plot_data <- long_spline %>% 
  filter(metabolite %in% metabolites_of_interest) %>% 
  left_join(names_and_class) %>% 
  bind_rows(long_crp %>%  mutate(metabolite = "CRP", class = "CRP") %>% select(metabolite, day, spline, class)) %>% 
  group_by(metabolite) %>% dplyr::mutate(y_label = spline[day == 30]) %>% ungroup() %>% 
  mutate(metabolite = ifelse(is.na(name), metabolite, name)) %>% 
  mutate(metabolite = ifelse(is.na(class), "PCT", metabolite)) %>% 
  mutate(class = ifelse(is.na(class), "PCT", class)) %>% 
  mutate(class = relevel(as.factor(class), "PCT"))

lpc_plot <- plot_data %>% 
  ggplot(aes(x = day, color = class)) +
  geom_line(aes(group = metabolite, y = spline, linetype = class), alpha = 0.7) +
  geom_text(aes(label = metabolite, x = 30, y = y_label), hjust = -0.02, show.legend = F, size = 2) +
  scale_linetype_manual(values = c(2, 3, 1, 1)) +
  scale_color_manual(values = c("black", "black", "#0563A9", "#71791D")) +
  scale_x_time_squish(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Time (days)", y = "Scaled metabolite and CRP levels", linetype = "", color = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2,byrow = TRUE), linetype = guide_legend(nrow = 2,byrow = TRUE))

#Gather figures
figure3b <- lpc_plot
tableS2 <- names_and_class %>% 
  filter(metabolite %in% attr(data.pretreated, "metrange")) %>%  
  select(-metabolite)
colnames(tableS2) <- c("Biochemical class", "Metabolite")

save(figure3b, tableS2, file = "manuscript/figures/plots_splines.Rdata")
