# Read results and plot mono-di mix ####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(grid)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)

mix_heptares <- readRDS("opp_simulationinfo_hepta_mixmonodi_7days.rds")
mix_heptares_b <- readRDS("opp_simulationinfo_hepta_mixmonodi_opp_7days.rds")

mix_W_hepta <- mix_heptares[[1]]
mix_res_hepta <- mix_heptares[[3]]
mix_W_hepta_est <- mix_res_hepta$exposures

mix_res_heptaopp <- mix_heptares_b[[3]]
mix_W_heptaopp_est <- mix_res_heptaopp$exposures




mix_data_hepta <- data.frame("True" = as.vector(mix_W_hepta), "Estimated" = as.vector(mix_W_hepta_est))
mix_p5 <- ggplot(data = mix_data_hepta, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("Hepta-nucleotide context - Mix mono-di")

mix_data_heptaopp <- data.frame("True" = as.vector(mix_W_hepta), "Estimated" = as.vector(mix_W_heptaopp_est))
mix_p6 <- ggplot(data = mix_data_heptaopp, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fill = "limegreen",
              color = "limegreen", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures")










# Read results and plot mono ####
#rm(list = ls())
library(ggplot2)
library(ggpubr)
library(grid)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
trires <- readRDS("opp_mono_simulationinfo_tri.rds")
W_tri <- trires[[1]]
res_tri <- trires[[3]]
W_tri_est <- res_tri$exposures

res_triopp <- trires[[4]]
W_triopp_est <- res_triopp$exposures

pentares <- readRDS("opp_mono_simulationinfo_penta.rds")
W_penta <- pentares[[1]]
res_penta <- pentares[[3]]
W_penta_est <- res_penta$exposures

res_pentaopp <- pentares[[4]]
W_pentaopp_est <- res_pentaopp$exposures

heptares <- readRDS("opp_mono_simulationinfo_heptaopp2.rds")
heptares_b <- readRDS("opp_mono_simulationinfo_heptaopp2_10init.rds")

W_hepta <- heptares[[1]]
res_hepta <- heptares_b[[3]]
W_hepta_est <- res_hepta$exposures

res_heptaopp <- heptares[[3]]
W_heptaopp_est <- res_heptaopp$exposures


data_tri <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_tri_est))
p1 <- ggplot(data = data_tri, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "#9D6ACF", method = "lm", se = TRUE, 
              color = "#9D6ACF", alpha = 0.3) + 
  ggtitle("Tri-nucleotide context") +
  ylab("Estimated exposures") +
  xlab("True exposures") + theme_minimal()

data_triopp <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_triopp_est))
p2 <- ggplot(data = data_triopp, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "#9D6ACF", 
              color = "#9D6ACF", alpha = 0.3)+
  ylab("Estimated exposures") +
  xlab("True exposures") + theme_minimal()


data_penta <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_penta_est))
p3 <- ggplot(data = data_penta, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "sandybrown", method = "lm", se = TRUE, 
              color = "sandybrown", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") + 
  ggtitle("Penta-nucleotide context") + theme_minimal()

data_pentaopp <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_pentaopp_est))
p4 <- ggplot(data = data_pentaopp, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "sandybrown", 
              color = "sandybrown", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") + theme_minimal()

data_hepta <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_hepta_est))
p5 <- ggplot(data = data_hepta, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("Hepta-nucleotide context") + theme_minimal()

data_heptaopp <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_heptaopp_est))
p6 <- ggplot(data = data_heptaopp, aes(x = True, y = Estimated)) +
  geom_point()  + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "palegreen2",
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") + theme_minimal()

void_plot <- textGrob(" ", 
                      gp = gpar(fontsize = 14, fontface = "bold"))
mid_title <- textGrob("Additive model with opportunities", 
                      gp = gpar(fontsize = 14, fontface = "bold"))

allplot0 <- ggarrange(p1,p3,p5,
                      void_plot, mid_title, void_plot,
                      p2,p4,p6,
                      nrow = 3, ncol = 3,
                      heights = c(1, 0.15, 1))
annotate_figure(allplot0, top = text_grob("Additive model without opportunities", face = "bold", size = 14))

allplot0 <- ggarrange(p1,p3,mix_p5,p5,
                      void_plot, mid_title, void_plot, void_plot,
                      p2,p4,mix_p6,p6,
                      nrow = 3, ncol = 4,
                      heights = c(1, 0.15, 1))
annotate_figure(allplot0, top = text_grob("Without opportunities", face = "bold", size = 14))


p5b <- ggplot(data = data_hepta, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("Without opportunities")

p6b <- ggplot(data = data_heptaopp, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fill = "limegreen",
              color = "limegreen", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("With opportunities")
allplot0 <- ggarrange(p5b,p6b)
annotate_figure(allplot0, top = text_grob("Hepta-nucleotide context", face = "bold", size = 14))

mhepta <- lm(Estimated ~0 +True, data_hepta)
summary(mhepta)
mheptaopp <- lm(Estimated ~0 + True, data_heptaopp, )
summary(mheptaopp)

mpentaopp <- lm(Estimated ~0 + True, data_pentaopp, )
summary(mpentaopp)

mtriopp <- lm(Estimated ~0 + True, data_triopp, )
summary(mtriopp)

mhepta <- lm(Estimated ~0 + True, data_hepta, )
summary(mhepta)

mpenta <- lm(Estimated ~0 + True, data_penta, )
summary(mpenta)

mtri <- lm(Estimated ~0 + True, data_tri, )
summary(mtri)
# Read results and plot di ####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(grid)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
trires <- readRDS("opp_simulationinfo_tri.rds")
W_tri <- trires[[1]]
res_tri <- trires[[3]]
W_tri_est <- res_tri$exposures

res_triopp <- trires[[4]]
W_triopp_est <- res_triopp$exposures

pentares <- readRDS("opp_simulationinfo_penta.rds")
W_penta <- pentares[[1]]
res_penta <- pentares[[3]]
W_penta_est <- res_penta$exposures

res_pentaopp <- pentares[[4]]
W_pentaopp_est <- res_pentaopp$exposures

heptares <- readRDS("opp_simulationinfo_heptaopp.rds")
W_hepta <- heptares[[1]]
res_hepta <- heptares[[3]]
W_hepta_est <- res_hepta$exposures

heptares <- readRDS("opp_simulationinfo_heptaopp_10init.rds")

res_heptaopp <- heptares[[3]]
W_heptaopp_est <- res_heptaopp$exposures


data_tri <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_tri_est))
p1 <- ggplot(data = data_tri, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "#9D6ACF", method = "lm", se = TRUE, 
              color = "#9D6ACF", alpha = 0.3) + 
  ggtitle("Tri-nucleotide context") +
  ylab("Estimated exposures") +
  xlab("True exposures") +theme_minimal()

data_triopp <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_triopp_est))
p2 <- ggplot(data = data_triopp, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "#9D6ACF", 
              color = "#9D6ACF", alpha = 0.3)+
  ylab("Estimated exposures") +
  xlab("True exposures")+theme_minimal()


mtri <- lm(Estimated ~0 +True, data_tri)
summary(mtri)
mtriopp <- lm(Estimated ~0 + True, data_triopp, )
summary(mtriopp) 

data_penta <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_penta_est))
p3 <- ggplot(data = data_penta, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "sandybrown", method = "lm", se = TRUE, 
              color = "sandybrown", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") + 
  ggtitle("Penta-nucleotide context") +theme_minimal()

data_pentaopp <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_pentaopp_est))
p4 <- ggplot(data = data_pentaopp, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "sandybrown", 
              color = "sandybrown", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures")+theme_minimal()


mpenta <- lm(Estimated ~True, data_penta)
summary(mpenta)
mpentaopp <- lm(Estimated ~True, data_pentaopp, )
summary(mpentaopp) 



data_hepta <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_hepta_est))
p5 <- ggplot(data = data_hepta, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("Hepta-nucleotide context")+theme_minimal()

data_heptaopp <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_heptaopp_est))
p6 <- ggplot(data = data_heptaopp, aes(x = True, y = Estimated)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = "gray", linewidth =1) +
  geom_smooth(method = "lm", se = TRUE, fill = "palegreen2",
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +theme_minimal()


void_plot <- textGrob(" ", 
                      gp = gpar(fontsize = 14, fontface = "bold"))
mid_title <- textGrob("Interaction model with opportunities", 
                      gp = gpar(fontsize = 14, fontface = "bold"))

allplot0 <- ggarrange(p1,p3,p5,
                      void_plot, mid_title, void_plot,
                      p2,p4,p6,
                      nrow = 3, ncol = 3,
                      heights = c(1, 0.15, 1))
annotate_figure(allplot0, top = text_grob("Interaction model without opportunities", face = "bold", size = 14))



mhepta <- lm(Estimated ~True, data_hepta)
summary(mhepta)
mheptaopp <- lm(Estimated ~ True, data_heptaopp, )
summary(mheptaopp) 
#add mono 
heptares_mono <- readRDS("opp_mono_simulationinfo_heptaopp2.rds")
heptares_b_mono <- readRDS("opp_mono_simulationinfo_heptaopp2_10init.rds")

W_hepta_mono <- heptares_mono[[1]]
res_hepta_mono <- heptares_b_mono[[3]]
W_hepta_est_mono <- res_hepta_mono$exposures

res_heptaopp_mono <- heptares_mono[[3]]
W_heptaopp_est_mono <- res_heptaopp_mono$exposures
data_hepta_mono <- data.frame("True" = as.vector(W_hepta_mono), "Estimated" = as.vector(W_hepta_est_mono))
p5_mono <- ggplot(data = data_hepta_mono, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") #+
#  ggtitle("Hepta-nucleotide context (Additive model)")

data_heptaopp_mono <- data.frame("True" = as.vector(W_hepta_mono), "Estimated" = as.vector(W_heptaopp_est_mono))
p6_mono <- ggplot(data = data_heptaopp_mono, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fill = "limegreen",
              color = "limegreen", alpha = 0.3) +
  ylab("") +
  xlab("True exposures")# +
  #ggtitle("Hepta-nucleotide context (Additive model)")


void_plot <- textGrob(" ", 
                     gp = gpar(fontsize = 14, fontface = "bold"))
# mid_title <- textGrob("With opportunities", 
#                       gp = gpar(fontsize = 14, fontface = "bold"))

library(cowplot)
allplot0 <- ggarrange(void_plot,p1,p3,p5,p5_mono, ncol = 1, nrow = 5,
                      heights = c(0.15, rep(1,4)))
tmp_p0<- ggdraw(allplot0) + 
  draw_label("Without opportunities", 
             x = 0.6, y = 0.975,      # Position (0 to 1)
             hjust = 0.5, vjust = 0.5, 
             fontface = 'bold')
allplot1 <- ggarrange(void_plot,p2,p4,p6,p6_mono, ncol = 1, nrow = 5,
                      heights = c(0.15, rep(1,4)))
tmp_p1 <- ggdraw(allplot1) + 
  draw_label("With opportunities", 
             x = 0.6, y = 0.975,      # Position (0 to 1)
             hjust = 0.5, vjust = 0.5, 
             fontface = 'bold')
ggarrange(tmp_p0, tmp_p1, ncol = 2, nrow = 1)

# Read results and plot full ####
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(grid)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
trires <- readRDS("opp_full_simulationinfo_tri.rds")
W_tri <- trires[[1]]
res_tri <- trires[[3]]
W_tri_est <- res_tri$exposures

res_triopp <- trires[[4]]
W_triopp_est <- res_triopp$exposures

pentares <- readRDS("opp_full_simulationinfo_penta.rds")
W_penta <- pentares[[1]]
res_penta <- pentares[[3]]
W_penta_est <- res_penta$exposures

res_pentaopp <- pentares[[4]]
W_pentaopp_est <- res_pentaopp$exposures

heptares <- readRDS("opp_full_simulationinfo_heptaopp3_10init.rds")
W_hepta <- heptares[[1]]
res_hepta <- heptares[[3]]
W_hepta_est <- res_hepta$exposures

#heptares2 <- readRDS("opp_full_simulationinfo_heptaopp.rds")
res_heptaopp <- heptares[[4]]
W_heptaopp_est <- res_heptaopp$exposures


data_tri <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_tri_est))
p1 <- ggplot(data = data_tri, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "orchid", method = "lm", se = TRUE, 
              color = "orchid", alpha = 0.3) + 
  ggtitle("Tri-nucleotide context") +
  ylab("Estimated exposures") +
  xlab("True exposures")

data_triopp <- data.frame("True" = as.vector(W_tri), "Estimated" = as.vector(W_triopp_est))
p2 <- ggplot(data = data_triopp, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fill = "mediumpurple3", 
              color = "mediumpurple3", alpha = 0.3)+
  ylab("Estimated exposures") +
  xlab("True exposures")


data_penta <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_penta_est))
p3 <- ggplot(data = data_penta, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "sandybrown", method = "lm", se = TRUE, 
              color = "sandybrown", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") + 
  ggtitle("Penta-nucleotide context") 

data_pentaopp <- data.frame("True" = as.vector(W_penta), "Estimated" = as.vector(W_pentaopp_est))
p4 <- ggplot(data = data_pentaopp, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fill = "darkorange1", 
              color = "darkorange1", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures")

data_hepta <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_hepta_est))
p5 <- ggplot(data = data_hepta, aes(x = True, y = Estimated)) +
  geom_point() +
  geom_smooth(fill = "palegreen2", method = "lm", se = TRUE,
              color = "palegreen2", alpha = 0.3) +
  ylab("Estimated exposures") +
  xlab("True exposures") +
  ggtitle("Hepta-nucleotide context")

 data_heptaopp <- data.frame("True" = as.vector(W_hepta), "Estimated" = as.vector(W_heptaopp_est))
 p6 <- ggplot(data = data_heptaopp, aes(x = True, y = Estimated)) +
   geom_point() +
   geom_smooth(method = "lm", se = TRUE, fill = "limegreen",
               color = "limegreen", alpha = 0.3) +
   ylab("Estimated exposures") +
   xlab("True exposures") #+
   # ylim(0,50000) +
   # xlim(0,25000)
   
void_plot <- textGrob(" ", 
                      gp = gpar(fontsize = 14, fontface = "bold"))
mid_title <- textGrob("With opportunities", 
                      gp = gpar(fontsize = 14, fontface = "bold"))

allplot0 <- ggarrange(p1,p3,p5,
                      void_plot, mid_title, void_plot,
                      p2,p4,p6,
                      nrow = 3, ncol = 3,
                      heights = c(1, 0.15, 1))
annotate_figure(allplot0, top = text_grob("Without opportunities", face = "bold", size = 14))

