#Put figures together
library(gridExtra)
library(tidyverse)


#Figures from correlations script and splines
load("manuscript/figures/plots_splines.Rdata")
load("manuscript/figures/plots_correlations.Rdata")

png(filename = "manuscript/figures/figure3.png", res = 600, width = 3.5*2.5, height = 5, units = "in")
grid.arrange(figure3a + labs(tag = "(A)"), figure3b + labs(tag = "(B)"), layout_matrix = t(c(1, 1, 1, 2, 2, 2, 2)))
dev.off()

png(filename = "manuscript/figures/figure5.png", res = 600, width = 10.05, height = 5, units = "in")
grid.arrange(figure5a + labs(tag = "(A)"), figure5b + labs(tag = "(B)"), layout_matrix = t(c(1, 2)))
dev.off()

png(filename = "manuscript/figures/figureS3.png", res = 300, width = 6, height = 9, units = "in")
print(figureS3)
dev.off()

png(filename = "manuscript/figures/figure4.png", res = 600, width = 6, height = 4, units = "in")
print(figure4)
dev.off()

write.csv(tableS2, file = "manuscript/tables/tableS2.csv")


#Figures from PCA script
load("manuscript/figures/plots_PCA.Rdata")

png(filename = "manuscript/figures/figure1.png", res = 600, width = 6, height = 2.6, units = "in")
print(figure1)
dev.off()


png(filename = "manuscript/figures/figureS1.png", res = 300, width = 0.75*9, height = 9, units = "in")
print(figureS1)
dev.off()

png(filename = "manuscript/figures/figureS2.png", res = 300, width = 0.75*9, height = 7, units = "in")
print(figureS2)
dev.off()

png(filename = "manuscript/figures/figure2.png", res = 600, width = 7, height = 11, units = "in")
grid.arrange(figure2a + labs(tag = "(A)"), figure2b + labs(tag = "(B)"), layout_matrix = matrix(c(1, 1, 2, 2, 2, 2)))
dev.off()

