#Exploratory analysis

# Load libraries
library(tidyverse)

#Metabolites
data.reduced <- read.csv2("data/00_data_reduced.csv")

#metabolites range
load("data/column_names.Rdata")

# Plot all metabolites
# of interest: subject.id, day
metabolite.data <- data.reduced[, c("subject.id", "day", metrange)]

metabolite.plot.df <- metabolite.data %>% 
  pivot_longer(-c(subject.id, day), names_to = "metabolite", values_to = "concentration") %>% 
  filter(!is.na(concentration))


squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(scales::trans_new("squished", trans, inv))
}


ticks <- 0:30
showing <- c(0:4,  30)
ticks[!ticks %in% showing] <- ""

all.metabolites <- unique(metabolite.plot.df$metabolite)
nr.of.pages <- 15
nr.of.metabolites <- 25
plot.metabolites <- list()
for (i in 1:nr.of.pages) {
  range.plots <- (1 + (i - 1) * nr.of.metabolites):(i * nr.of.metabolites)
  range.plots <- range.plots[range.plots <= length(all.metabolites)]
  
  plot.metabolites[[i]] <- metabolite.plot.df %>%
    filter(metabolite %in% all.metabolites[range.plots]) %>% 
    ggplot(aes(x = day, y = concentration, color = as.factor(subject.id), group = as.factor(subject.id))) +
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) + 
    scale_x_continuous(trans = squish_trans(4.5, 29.5, 7), breaks = 0:30, labels = as.character(ticks)) +
    
    facet_wrap(~ metabolite, scale = "free_y") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
} 


pdf(file = "results/figures/metabolomics_over_time.pdf", width = 11, height = 11)
for (i in 1:length(plot.metabolites)) {
  print(plot.metabolites[[i]])
}
dev.off()

  