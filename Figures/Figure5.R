rm(list=ls())

## BIC plot 

# load output results
#load("BRCAres3_5ThreeInitNBwithOPP.RData")
load("BRCA_3_5ThreeInitNBwithOPP.RData")
brcares35 = brcares

#load("BRCAres7ThreeInitNBwithOPP.RData")
load("BRCA_7ThreeInitNBwithOPP.RData")
brcares <- c(brcares35[c(1,2,3)], brcares[2], brcares35[c(4,5)], brcares[3],
             brcares35[c(6,7)], brcares[4])

names(brcares)[length(names(brcares))] <- "full7opp"

data_opp = brcares

# load output results
#load("BRCAres3_5ThreeInitNBwithoutOPP.RData")
load("BRCA_3_5ThreeInitNBwithoutOPP.RData")
brcares35 = brcares

#load("BRCAres7ThreeInitNBwithoutOPP.RData")
load("BRCA_7ThreeInitNBwithoutOPP.RData")
brcares <- c(brcares35[c(1,2,3)], brcares[2], brcares35[c(4,5)], brcares[3],
             brcares35[c(6,7)], brcares[4])

names(brcares)[length(names(brcares))] <- "full7"

data = brcares

varnames = names(brcares)

load("data/BRCA/brcacounts3.RData")
load("data/BRCA/brcacounts5.RData")
load("data/BRCA/brcacounts7.RData")
load("data/opportunities.RData")
navn = c("Additive","Interaction","Standard")
navn_no = c("Tri","Penta","Hepta")
no = c(3,5,7)
param = c("mono","didi","full")
signatures = 4

library(ggplot2)
data_all = c()
opp_list = list(opp3,opp5,opp7)
data_list = list(brcacounts3,brcacounts5,brcacounts7)

plot_list = list()
for(i in c(1:3)){
  for(j in c(1:3)){
    opp = opp_list[[j]]
    
    estimate_opp = data_opp[[paste0(param[i],no[j],"opp")]][[signatures/2]]$exposures%*%data_opp[[paste0(param[i],no[j],"opp")]][[signatures/2]]$signatures
    estimate_opp = t(t(estimate_opp)*opp)
    estimate = data[[paste0(param[i],no[j])]][[signatures/2]]$exposures%*%data[[paste0(param[i],no[j])]][[signatures/2]]$signatures

  
  res_opp = ((estimate_opp + 1e-10) - data_list[[j]])^2/(estimate_opp + 1e-10)
  res = ((estimate + 1e-10) - data_list[[j]])^2/(estimate + 1e-10)
  #res_opp = estimate_opp - data_list[[j]]*log(estimate_opp + 1e-10)
  #res = estimate - data_list[[j]]*log(estimate + 1e-10)

  diff =  (colMeans(res) - colMeans(res_opp))

  id = order(opp)

  df = data.frame(id = 1:length(opp),opp = opp[id], diff = diff[id])
  df = df[df$id < quantile(df$id, 0.1) | df$id > quantile(df$id, 0.9),]
  divide = 10
  repeach = floor(nrow(df)/divide)
  left = nrow(df) %% repeach
  df$groups = c(rep(c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'), each = repeach),rep('90-100',left))
  #df$groups = c(rep(1:divide, each = repeach),rep(divide,left))
  print(df$groups[1:10])
  df$groups = as.factor(df$groups)
  
  df$param = navn[i]
  df$context = navn_no[j]
  data_all = rbind(data_all,df)
  # maxq = quantile(df$diff, 0.99)
  # df$diff[df$diff > maxq] = maxq
  # minq = quantile(df$diff, 0.01)
  # df$diff[df$diff < minq] = minq
  if(i == 1){
  plot_list[[paste0(param[i],no[j])]] = ggplot(df, aes(x = groups,y = diff))+
    #geom_bar(stat="identity")+
    geom_hline(yintercept = 0, col = "red", lwd = 1)+
    geom_boxplot(outlier.shape = NA)+
    # geom_jitter(cex = 1)+
    
    # geom_smooth(method = "loess")+
    # ylim(-quantile(df$diff,0.93),quantile(df$diff,0.93))+
    # scale_fill_gradient2(low="blue",mid = "white", high="red", midpoint = 0, name = "prediction improvement")+
    #guides(fill = guide_legend(title = "log prediction improvement"))+
    theme_bw()+
    theme(legend.title = element_text(size = 5),
          #axis.text = element_blank(), 
          #axis.ticks = element_blank(), 
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size=5),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 90,hjust=0.95, vjust = -0.1)
          #panel.background = element_blank()
          )+
    # scale_x_log10()+
    ylab(paste0("\n Standardized difference"))+
    
    xlab("Opportunity")+
    ggtitle(navn_no[j])
  }else{
  plot_list[[paste0(param[i],no[j])]] = ggplot(df, aes(x = groups,y = diff))+
    #geom_bar(stat="identity")+
    geom_hline(yintercept = 0, col = "red", lwd = 1)+
    geom_boxplot(outlier.shape = NA)+
    # geom_jitter(cex = 1)+
    # geom_point(cex = 1)+
    
    # geom_smooth(method = "loess")+
    # ylim(-quantile(df$diff,0.95),quantile(df$diff,0.95))+
    #scale_fill_gradient2(low="blue",mid = "white", high="red", midpoint = 0, name = "log prediction improvement")+
    #guides(fill = guide_legend(title = "log prediction improvement"))+
    theme_bw()+
    theme(legend.title = element_text(size = 5),
          #axis.text = element_blank(), 
          #axis.ticks = element_blank(), 
          legend.key.size = unit(0.5, 'cm'),
          legend.text = element_text(size=5),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 90,hjust=0.05,vjust=0.2)
          #panel.background = element_blank()
      )+
      # scale_x_log10()+
      ylab(paste0("\n Standardized difference"))+
      xlab("Opportunity")+
      ggtitle("")
}
}
}


library(ggpubr)

p1 = ggarrange(plotlist = plot_list, nrow = 3,ncol = 3)
p1
# annotate_figure(p1, top = text_grob("Hepta", face = "bold", size = 14))
annotate_figure(p1, left = text_grob("Standard                      Interaction                      Additive", rot = 90),)

data_all$context = factor(data_all$context, levels = c("Tri","Penta","Hepta"))

fills = c('#AC0136',"#FF9912","#27408B")
data_all$color_y = '#AC0136'
data_all$color_y[data_all$param == 'Interaction'] = "#FF9912"
data_all$color_y[data_all$param == 'Standard'] = "#27408B"
ggplot(data_all, aes(x = groups, y = diff))+
  geom_hline(yintercept = 0, col = "red", lwd = 1)+
  geom_boxplot(width = 0.8,outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  #geom_boxplot(outlier.shape = NA)+
  #gg.layers::geom_boxplot2(width = 0.8, width.errorbar = 0.5)+
  # geom_blank(data=dummy) +
  ylim(-3,3)+
  #ylim(-quantile(df$diff,0.95),quantile(df$diff,0.95))+
  #facet_grid(param ~ context, scales = "free_y")+
  ggh4x::facet_grid2(param ~ context, scales = "free_y", independent = "y")+
  theme_bw()+
  theme(legend.title = element_text(size = 5),
        #strip.background.y = element_rect(fill=fills),
        # strip.text.y = element_text(color = data_all$color_y),
        #axis.text = element_blank(), 
        #axis.ticks = element_blank(), 
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=5),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5)
        #panel.background = element_blank()
  )+
  # scale_x_log10()+
  ylab(paste0("Residual difference"))+
  xlab("Opportunity percentiles")

head(data_all)
