# fig 2 ####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(grid)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)

setwd("/home/au687642/OneDrive/Documents/Postdoc/MutationalOpportunities/ParOpp/data/")
load("opportunities.RData")
load("BRCA/brcacounts3.RData")
load("BRCA/brcacounts5.RData")
load("BRCA/brcacounts7.RData")
load("LIRI/livercounts3.RData")
load("LIRI/livercounts5.RData")
load("LIRI/livercounts7.RData")



p1_data <- data.frame("Counts" = colMeans(brcacounts3), "Opportunities" = opp3)
p1 <- ggplot(data = p1_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "#9D6ACF", method = "lm", se = TRUE,
              color = "#9D6ACF", alpha = 0.3) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts") +
  ggtitle("Tri-nucleotide context")

p2_data <- data.frame("Counts" = colMeans(brcacounts5), "Opportunities" = opp5)
p2 <- ggplot(data = p2_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "sandybrown", method = "lm", se = TRUE,
              color = "sandybrown", alpha = 0.3) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts") +
  ggtitle("Penta-nucleotide context")

p3_data <- data.frame("Counts" = colMeans(brcacounts7), "Opportunities" = opp7)
p3 <- ggplot(data = p3_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "limegreen", method = "lm", se = TRUE,
              color = "limegreen", alpha = 0.3) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts") +
  ggtitle("Hepta-nucleotide context")



p4_data <- data.frame("Counts" = colMeans(livercounts3), "Opportunities" = opp3)
p4 <- ggplot(data = p4_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "#9D6ACF", method = "lm", se = TRUE,
              color = "#9D6ACF", alpha = 0.3)+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts")

p5_data <- data.frame("Counts" = colMeans(livercounts5), "Opportunities" = opp5)
p5 <- ggplot(data = p5_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "sandybrown", method = "lm", se = TRUE,
              color = "sandybrown", alpha = 0.3)+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts") 

p6_data <- data.frame("Counts" = colMeans(livercounts7), "Opportunities" = opp7)
p6 <- ggplot(data = p6_data, aes(x = Counts, y = Opportunities)) +
  geom_point() +
  geom_smooth(fill = "limegreen", method = "lm", se = TRUE,
              color = "limegreen", alpha = 0.3) +
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0) +
  xlab("Mutational counts") 

allplot <- ggarrange(p1,p2,p3,
                      p4,p5,p6,
                      nrow = 2, ncol = 3)
allplot
