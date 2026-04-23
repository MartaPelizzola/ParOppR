# Read results BRCA 3-5 OPP with percentile ####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(SigMoS)
load("data/BRCA/brcacounts3.RData")
load("data/BRCA/brcacounts5.RData")
load("data/BRCA/brcacounts7.RData")

load("data/opportunities.RData")
load("BRCA_3_5ThreeInitNBwithOPP.RData")
brca3_5_opp <- brcares
load("BRCA_3_5ThreeInitNBwithoutOPP.RData")
# Three without opp ####
est <- (brcares$didi3[[2]]$exposures%*%brcares$didi3[[2]]$signatures)
rsd_nb3 = (brcacounts3 - est)/sqrt(est+1e-10)

dataNB = data.frame(obs_nb = opp3, Residuals_nb = colMeans(rsd_nb3))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r1 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() + 
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x,#bs(x, df = 2), 
                color = "orchid", 
                size = 1) +
  scale_x_continuous(breaks=(0:9)/9,
                   labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  #scale_x_log10()+
  #ylim(-500,500) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         
        axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("") +
  ggtitle("Without opportunities")

# Three with opp ####
est_opp <- (brca3_5_opp$didi3opp[[4]]$exposures%*%brca3_5_opp$didi3opp[[4]]$signatures)*matrix(rep(opp3, nrow(brcacounts3)), ncol = length(opp3), nrow = nrow(brcacounts3), byrow = T)
rsd_nb3_opp = (brcacounts3 - est_opp)/sqrt(est_opp+1e-10)

dataNB = data.frame(obs_nb = opp3, Residuals_nb = colMeans(rsd_nb3_opp))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r2 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() +
  ylim(-1,1) +
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x,#bs(x, df = 2), 
                color = "mediumpurple3", 
                size = 1) +
  #geom_quantile(quantiles = 0.95, color = "darkgreen", linewidth = 1) +
  #geom_quantile(quantiles = 0.05, color = "darkgreen", linewidth = 1) +
  #scale_x_log10()+
  #  ylim(-500,500) +
  scale_x_continuous(breaks=(0:9)/9,
                     labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),         
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("")+
  ggtitle("With opportunities")

p_three0 <- ggarrange(r1,r2,ncol = 2)
p_three <- annotate_figure(p_three0, right = text_grob("Tri-nucleotide", rot = 90, size = 20, hjust = 0.2))

# Five without opp ####
est <- (brcares$didi5[[4]]$exposures%*%brcares$didi5[[4]]$signatures)
rsd_nb5 = (brcacounts5 - est)/sqrt(est+1e-10)

dataNB = data.frame(obs_nb = opp5, Residuals_nb = colMeans(rsd_nb5))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r3 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() + 
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x,#bs(x, df = 2), 
                color = "sandybrown", 
                size = 1) +
  #scale_x_log10()+
  ylim(-1,3) +
  scale_x_continuous(breaks=(0:9)/9,
                     labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         
        axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),         
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("") 
# Five with opp ####
est_opp <- (brca3_5_opp$didi5opp[[4]]$exposures%*%brca3_5_opp$didi5opp[[4]]$signatures)*matrix(rep(opp5, nrow(brcacounts5)), ncol = length(opp5), nrow = nrow(brcacounts5), byrow = T)
rsd_nb5_opp = (brcacounts5 - est_opp)/sqrt(est_opp+1e-10)

dataNB = data.frame(obs_nb = opp5, Residuals_nb = colMeans(rsd_nb5_opp))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r4 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() +
  ylim(-1,3) +
  scale_x_continuous(breaks=(0:9)/9,
                     labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x,#bs(x, df = 2), 
                color = "darkorange1", 
                size = 1) +
  #geom_quantile(quantiles = 0.95, color = "darkgreen", linewidth = 1) +
  #geom_quantile(quantiles = 0.05, color = "darkgreen", linewidth = 1) +
  #scale_x_log10()+
#  ylim(-500,500) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         
        axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),         
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("")


p_five0 <- ggarrange(r3,r4,ncol = 2)
p_five <- annotate_figure(p_five0, right = text_grob("Penta-nucleotide", rot = 90, size = 20, hjust = 0.2))

# Seven without opp ####
load("BRCA_7ThreeInitNBwithoutOPP.RData")
est <- (brcares$didi7[[4]]$exposures%*%brcares$didi7[[4]]$signatures)
rsd_nb7 = (brcacounts7 - est)/sqrt(est+1e-10)

dataNB = data.frame(obs_nb = opp7, Residuals_nb = colMeans(rsd_nb7))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r5 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() + 
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x, #bs(x, df = 2), 
                color = "palegreen2", 
                size = 1) +
  #scale_x_log10()+
  ylim(-1,5) +
  scale_x_continuous(breaks=(0:9)/9,
                     labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         
        axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),         
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("") 
# Seven with opp ####
load("BRCA_7ThreeInitNBwithOPP.RData")
est_opp <- (brcares$didi7opp[[4]]$exposures%*%brcares$didi7opp[[4]]$signatures)*matrix(rep(opp7, nrow(brcacounts7)), ncol = length(opp7), nrow = nrow(brcacounts7), byrow = T)
rsd_nb7_opp = (brcacounts7 - est_opp)/sqrt(est_opp+1e-10)

dataNB = data.frame(obs_nb = opp7, Residuals_nb = colMeans(rsd_nb7_opp))
dataNB <- dataNB %>%
  mutate(opp_percentile = percent_rank(obs_nb))

r6 <- ggplot(dataNB, aes(x = opp_percentile, y = Residuals_nb))+
  geom_point() +
  geom_quantile(quantiles = c(0.05, 0.95), 
                formula = y ~ x,#bs(x, df = 2), 
                color = "limegreen", 
                size = 1) +
  #geom_quantile(quantiles = 0.95, color = "darkgreen", linewidth = 1) +
  #geom_quantile(quantiles = 0.05, color = "darkgreen", linewidth = 1) +
  #scale_x_log10()+
  ylim(-1,5) +
  scale_x_continuous(breaks=(0:9)/9,
                     labels=c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),         
        axis.text.x = element_text(angle = 90,hjust=0.95,size=14, vjust = 0.5),         
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(ggplot2::margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("") + ylab("")


p_Seven0 <- ggarrange(r5,r6,ncol = 2)
p_Seven <- annotate_figure(p_Seven0, right = text_grob("Hepta-nucleotide", rot = 90, size = 20, hjust = 0.2))

p_all0 <- ggarrange(p_three,p_five,p_Seven, ncol = 1, nrow = 3)
p_all <- annotate_figure(p_all0, bottom = text_grob("Opportunity percentile", size = 20, vjust = -1.5))
annotate_figure(p_all, left = text_grob("Standardized residuals", size = 20, rot = 90, vjust = 2))

