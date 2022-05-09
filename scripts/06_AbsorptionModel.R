#create absorption model for all metabolites with individual metabolite/patient parameters
library(tidyverse)

data.raw <- read.csv2("data/00_data_reduced.csv")
load("data/column_names.Rdata")
#data preparation to create NONMEM format
#columns: ID, TIME, DV, AMT, EVID, CMT, WT, OCC

long_data <- data.raw %>% pivot_longer(cols = all_of(metrange), names_to = "metabolite")# %>% 
#  filter(metabolite %in% metrange[1:3]) #filter only few for testing

data.frame(occasion_value = as.numeric(as.factor(metrange)), metabolite = metrange) %>% View

nonmem_data <- data.frame(ID = long_data$subject.id,
                          TIME = long_data$day,
                          DV = long_data$value,
                          AMT = 0,
                          CMT = 2,
                          EVID = 0,
                          OCC = as.numeric(as.factor(long_data$metabolite)),
                          TID = as.numeric(as.factor(paste0(long_data$subject.id, as.factor(long_data$metabolite))))) %>% 
  arrange(ID, OCC, TIME)


nonmem_rows <- nonmem_data %>% select(-c(TIME, DV)) %>% unique() %>% 
  mutate(TIME = 0, CMT = 1, DV = NA, EVID = 4, AMT = 1) %>% 
  select(names(nonmem_data))

nonmem_data_final <- rbind(nonmem_data, nonmem_rows) %>% as.data.frame() %>% 
  arrange(TID, TIME, -EVID) %>% 
  mutate(EVID = ifelse(AMT == 0 & EVID == 4,  3, EVID))

#one interesting metabolite: TG.52.3 - OCC 299
#exploratory plot
nonmem_data_final %>% 
  filter(OCC == 299) %>% 
  ggplot(aes(x = TIME, y = DV, group = ID, color = as.factor(ID))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ OCC, scales = "free_y") +
  theme_bw()

#absorption model
absortion_function <- function(t, intercept, ka, ke) {
  intercept*ka/(ka - ke)*(exp(-ke*t) - exp(-ka*t))
}
#derivative of absorption model for fitting nlmer model
absorption_derivative <- deriv(body(absortion_function)[[2]], 
                              namevec = c("intercept", "ka", "ke"), 
                              function.arg = absortion_function)

#filter the TG.52.3 metabolite
test_data <- nonmem_data_final %>% filter(OCC == 299)

#fit absorption model
fit <- lme4::nlmer(DV ~ absorption_derivative(t = TIME,  intercept, ka, ke) ~
                     (intercept | ID) + (ka | ID) + (ke | ID), data = test_data, 
                   start = c(intercept = 40, ka = 2.4, ke = 0.1))

#extract rendaom and fixed effects
ranef(fit)$ID + data.frame(t(fixef(fit)))[col(ranef(fit)$ID)]

#predict curve (fixed effect due to zero random effect values)
days <- seq(0, 30, by = 0.1)
pred_dat <- data.frame(TIME = days, DV = absortion_function(days, fixef(fit)[1], fixef(fit)[2], fixef(fit)[3]))

#plot both raw data and fixed model fit
pdf("results/figures/absorption_model_fit_TG.52.3.pdf", height = 7, width = 10)
nonmem_data_final %>% 
  filter(OCC == 299) %>% 
  ggplot(aes(x = TIME, y = DV, group = ID, color = as.factor(ID))) +
  geom_point() +
  geom_line() +
  geom_line(data = pred_dat, aes(x = TIME, y = DV), inherit.aes = FALSE) +
  facet_wrap(~ ID) +
  theme_bw()
dev.off()




