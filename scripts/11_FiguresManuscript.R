#Put figures together
library(gridExtra)
library(tidyverse)


#Figures from Kmeans script
load("manuscript/figures/plots_Kmeans.Rdata")

png(filename = "manuscript/figures/figure2.png", res = 600, width = 6, height = 3.5, units = "in")
print(figure2)
dev.off()

png(filename = "manuscript/figures/figureS1.png", res = 300, width = 7, height = 4.5, units = "in")
print(figureS1)
dev.off()

png(filename = "manuscript/figures/figureS2.png", res = 300, width = 7, height = 9, units = "in")
print(figureS2)
dev.off()

png(filename = "manuscript/figures/figureS3.png", res = 300, width = 10, height = 8.5, units = "in")
print(figureS3)
dev.off()

write.csv(tableS2, file = "manuscript/tables/tableS2.csv")

#Figures from correlations script
load("manuscript/figures/plots_correlations.Rdata")

png(filename = "manuscript/figures/figure3.png", res = 600, width = 3.5*2.5, height = 5, units = "in")
grid.arrange(figure3a + labs(tag = "(A)"), figure3b + labs(tag = "(B)"), layout_matrix = t(c(1, 1, 1, 2, 2, 2, 2)))
dev.off()

png(filename = "manuscript/figures/figure4.png", res = 600, width = 10.05, height = 5, units = "in")
grid.arrange(figure4a + labs(tag = "(A)"), figure4b + labs(tag = "(B)"), layout_matrix = t(c(1, 2)))
dev.off()

png(filename = "manuscript/figures/figureS4.png", res = 300, width = 6, height = 9, units = "in")
print(figureS4)
dev.off()

png(filename = "manuscript/figures/figureS5.png", res = 600, width = 4.3, height = 5.5, units = "in")
print(figureS5)
dev.off()


#Figures from PCA script
load("manuscript/figures/plots_PCA.Rdata")

png(filename = "manuscript/figures/figure6.png", res = 600, width = 6, height = 6, units = "in")
grid.arrange(figure6a + labs(tag = "A"), figure6b + labs(tag = "B"), figure6c+ labs(tag = "C"), layout_matrix = rbind(c(1,1), c(2,3)))
dev.off()

#
png(filename = "manuscript/figures/figure6a.png", res = 600, width = 6, height = 3, units = "in")
print(figure6a)
dev.off()


png(filename = "manuscript/figures/figureS6.png", res = 300, width = 0.75*9, height = 9, units = "in")
print(figureS6)
dev.off()





