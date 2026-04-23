rm(list = ls())
setwd("~/OneDrive/Documents/Postdoc/MutationalOpportunities/ParOpp/results/negbin/")
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)
library(dplyr)
############################ Functions ###############
matchSig <- function(H1, H2 = NULL){
  if (is.null(H2)){
    H2 <- read.delim("http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", header = TRUE)
    rownames(H2) <- H2$Somatic.Mutation.Type
    H2 <- t(H2[,grep("Signature", colnames(H2))])
  } else if (ncol(H1)!=ncol(H2)) {
    stop("The two signatures matrices must have the same number of mutation types.")
  }


  dist <- cor(t(H1),t(H2))
  K <- nrow(H1)
  m1 <- numeric(K)
  m2 <- numeric(K)
  distmat = dist
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- -2
    dist[,remove[1,2]] <- -2
    m1[s] <- remove[1,1]
    m2[s] <- remove[1,2]
  }
  Output <- list()
  Output$S1 <- H1[m1,]  # the best matched signatures from H1
  Output$S2 <- H2[m2,]  # the corresponding best matched signatures from H2
  Output$distmat <- distmat      # the distance matrix of correlations

  return(Output)
}
gkl.dev = function(y, mu){
  r = mu
  p = which(y > 0)
  r[p] = (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}

NBlik = function(M,P,E, alpha){
  y = as.vector(M)
  mu = as.vector(P%*%E)
  r = lgamma(y + alpha) - lfactorial(y) - lgamma(alpha) + alpha * log(alpha/(alpha + mu))
  p = which(y > 0)
  r[p] = (lgamma(y + alpha) - lfactorial(y) - lgamma(alpha) + y * log(mu/(alpha + mu)) + alpha * log(alpha/(alpha + mu)))[p]

  obj = sum(r) # euclidean distance

  return(obj)
}


# load data ####

load("../../data/BRCA/brcacounts3.RData")
load("../../data/BRCA/brcacounts5.RData")
load("../../data/BRCA/brcacounts7.RData")
ordermut3 <- order(substr(colnames(brcacounts3),3,5))
brcacounts3 <- brcacounts3[,order(substr(colnames(brcacounts3),3,5))]
brcacounts5 <- brcacounts5[,order(substr(colnames(brcacounts5),4,6))]
brcacounts7 <- brcacounts7[,order(substr(colnames(brcacounts7),5,7))]

cosmic <- read_table("~/OneDrive/Documents/Postdoc/MutationalScanning/Trios/COSMIC_v3.3.1_SBS_GRCh37.txt")
cosmic <- as.data.frame(cosmic)
rownames(cosmic) <- cosmic$Type
cosmic <- t(cosmic[,grep("SBS", colnames(cosmic))])
cosmic <- cosmic[,colnames(brcacounts3)]
cosmic <- cosmic[1:40,]
### load results without opp#########

load("BRCAres3_5ThreeInitNBwithoutOPP.RData")
brcares35 = brcares

load("BRCAres7ThreeInitNBwithoutOPP.RData")
brcares <- c(brcares35[c(1,2,3)], brcares[2], brcares35[c(4,5)], brcares[3],
             brcares35[c(6,7)], brcares[4])
names(brcares)[length(names(brcares))] <- "full7"
varnames <- names(brcares)

# cosine similarity ####
cossim <- list()
bestsig <- c()
for(i in varnames[c(2,5,8)]){
  rownames(brcares[[i]][[2]]$signatures) <- paste0("S", 1:4)
  brcares[[i]][[2]]$signatures <- brcares[[i]][[2]]$signatures[, ordermut3]
  cossim[[i]] <- matchSig(brcares[[i]][[2]]$signatures, cosmic)
  tmp <- matchSig(brcares[["mono3"]][[2]]$signatures, brcares[[i]][[2]]$signatures)
  bestsig <- cbind(bestsig, rownames(cossim[[i]]$S2)[match(rownames(tmp$S1),rownames(cossim[[i]]$S1))])
}


# plot ####
H <- list()
H_named <- list()
p <- list()
idx <- 1
sig <- 4
modelname <- c("Additive", "Interaction", "Standard")
colorsbs <- palette("ggplot2")[1:length(unique(c(rownames(cossim[[varnames[2]]]$S2),
                                                 rownames(cossim[[varnames[5]]]$S2),
                                                 rownames(cossim[[varnames[8]]]$S2))))]
names(colorsbs) <- unique(c(rownames(cossim[[varnames[2]]]$S2),rownames(cossim[[varnames[5]]]$S2),
                            rownames(cossim[[varnames[8]]]$S2)))
for (i in varnames[c(2,5,8)]){
  H[[i]] <- data.frame(t(cossim[[i]]$S1))
  colnames(H[[i]]) <- paste0(rownames(cossim[[i]]$S2), rep(" (",4), round(cossim[[i]]$distmat[cbind(rownames(cossim[[i]]$S1), rownames(cossim[[i]]$S2))],2), ")")
  rownames(H[[i]]) <- colnames(brcacounts3)
  H_named[[i]] <- H[[i]]
  H[[i]]$MutationType <- colnames(brcacounts3)
  H[[i]]$Substitution = factor(str_sub(H[[i]]$MutationType, 3, -3), levels = unique(str_sub(H[[i]]$MutationType, 3, -3)))
  H[[i]] = pivot_longer(H[[i]], colnames(H[[i]])[1]:colnames(H[[i]])[sig])
  colnames(H[[i]])[c(3,4)] = c("Signature", "Intensity")
  H[[i]] = H[[i]] %>% arrange(Signature,Substitution)
  H[["mono3"]]$Signature <- factor(H[["mono3"]]$Signature,
                                   levels = c("SBS6 (0.87)", "SBS2 (0.81)",
                                              "SBS3 (0.57)","SBS25 (0.27)"))
  H[["didi3"]]$Signature <- factor(H[["didi3"]]$Signature,
                                   levels = c("SBS1 (0.9)", "SBS2 (0.78)",
                                              "SBS5 (0.75)", "SBS13 (0.65)"))
  H[["full3"]]$Signature <- factor(H[["full3"]]$Signature,
                                   levels = c("SBS1 (0.82)", "SBS2 (0.99)",
                                              "SBS3 (0.82)", "SBS13 (0.82)"))
  H[[i]]$MutationType <- factor(H[[i]]$MutationType, levels = unique(H[[i]]$MutationType))
  p[[idx]] <- ggplot(H[[i]]) + geom_col(aes(x = MutationType, y = Intensity, fill = Substitution))  + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 40),
          axis.title = element_text(size = 30),
          axis.text.y = element_text(size = 20),
          strip.text = element_text(size = 30),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30),
          strip.background = element_rect(colour="white", fill="white") ) + ylim(0,0.4) +
    facet_wrap(vars(Signature), nrow = 4, strip.position = "top",
               scales="free_x") + ggtitle(modelname[idx])
  idx = idx + 1
}

ggarrange(p[[1]],p[[2]],p[[3]], ncol = 3, common.legend = T, legend = "bottom")



##############################################################################
# OPPORTUNITIES ####
##############################################################################
rm(list = ls())
setwd("~/OneDrive/Documents/Postdoc/MutationalOpportunities/ParOpp/results/negbin/")
library(ggplot2)
library(ggpubr)
library(readr)
library(stringr)
library(tidyverse)
############################ Functions ###############
matchSig <- function(H1, H2 = NULL){
  if (is.null(H2)){
    H2 <- read.delim("http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", header = TRUE)
    rownames(H2) <- H2$Somatic.Mutation.Type
    H2 <- t(H2[,grep("Signature", colnames(H2))])
  } else if (ncol(H1)!=ncol(H2)) {
    stop("The two signatures matrices must have the same number of mutation types.")
  }
  
  
  dist <- cor(t(H1),t(H2))
  K <- nrow(H1)
  m1 <- numeric(K)
  m2 <- numeric(K)
  distmat = dist
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- -2
    dist[,remove[1,2]] <- -2
    m1[s] <- remove[1,1]
    m2[s] <- remove[1,2]
  }
  Output <- list()
  Output$S1 <- H1[m1,]  # the best matched signatures from H1
  Output$S2 <- H2[m2,]  # the corresponding best matched signatures from H2
  Output$distmat <- distmat      # the distance matrix of correlations
  
  return(Output)
}
gkl.dev = function(y, mu){
  r = mu
  p = which(y > 0)
  r[p] = (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}

NBlik = function(M,P,E, alpha){
  y = as.vector(M)
  mu = as.vector(P%*%E)
  r = lgamma(y + alpha) - lfactorial(y) - lgamma(alpha) + alpha * log(alpha/(alpha + mu))
  p = which(y > 0)
  r[p] = (lgamma(y + alpha) - lfactorial(y) - lgamma(alpha) + y * log(mu/(alpha + mu)) + alpha * log(alpha/(alpha + mu)))[p]
  
  obj = sum(r) # euclidean distance
  
  return(obj)
}


# load data ####

load("../../data/BRCA/brcacounts3.RData")
load("../../data/BRCA/brcacounts5.RData")
load("../../data/BRCA/brcacounts7.RData")
ordermut3 <- order(substr(colnames(brcacounts3),3,5))
brcacounts3 <- brcacounts3[,order(substr(colnames(brcacounts3),3,5))]
brcacounts5 <- brcacounts5[,order(substr(colnames(brcacounts5),4,6))]
brcacounts7 <- brcacounts7[,order(substr(colnames(brcacounts7),5,7))]

cosmic <- read_table("~/OneDrive/Documents/Postdoc/MutationalScanning/Trios/COSMIC_v3.3.1_SBS_GRCh37.txt")
cosmic <- as.data.frame(cosmic)
rownames(cosmic) <- cosmic$Type
cosmic <- t(cosmic[,grep("SBS", colnames(cosmic))])
cosmic <- cosmic[,colnames(brcacounts3)]
cosmic <- cosmic[1:40,]
### load results with opp#########

load("BRCAres3_5ThreeInitNBwithOPP.RData")
brcares35 = brcares

load("BRCAres7ThreeInitNBwithOPP.RData")
brcares <- c(brcares35[c(1,2,3)], brcares[2], brcares35[c(4,5)], brcares[3],
             brcares35[c(6,7)], brcares[4])

varnames = names(brcares)
# cosine similarity ####
cossim <- list()
bestsig <- c()
for(i in varnames[c(2,5,8)]){
  rownames(brcares[[i]][[2]]$signatures) <- paste0("S", 1:4)
  brcares[[i]][[2]]$signatures <- brcares[[i]][[2]]$signatures[, ordermut3]
  cossim[[i]] <- matchSig(brcares[[i]][[2]]$signatures, cosmic)
  tmp <- matchSig(brcares[["mono3opp"]][[2]]$signatures, brcares[[i]][[2]]$signatures)
  bestsig <- cbind(bestsig, rownames(cossim[[i]]$S2)[match(rownames(tmp$S1),rownames(cossim[[i]]$S1))])
}


# plot ####
H <- list()
H_named <- list()
p <- list()
idx <- 1
sig <- 4
modelname <- c("Additive", "Interaction", "Standard")
colorsbs <- palette("ggplot2")[1:length(unique(c(rownames(cossim[[varnames[2]]]$S2),
                                                 rownames(cossim[[varnames[5]]]$S2), 
                                                 rownames(cossim[[varnames[8]]]$S2))))]
names(colorsbs) <- unique(c(rownames(cossim[[varnames[2]]]$S2),rownames(cossim[[varnames[5]]]$S2),
                            rownames(cossim[[varnames[8]]]$S2))) 
for (i in varnames[c(2,5,8)]){
  H[[i]] <- data.frame(t(cossim[[i]]$S1))
  colnames(H[[i]]) <- paste0(rownames(cossim[[i]]$S2), rep(" (",4), round(cossim[[i]]$distmat[cbind(rownames(cossim[[i]]$S1), rownames(cossim[[i]]$S2))],2), ")")
  rownames(H[[i]]) <- colnames(brcacounts3)
  H_named[[i]] <- H[[i]]
  H[[i]]$MutationType <- colnames(brcacounts3)
  H[[i]]$Substitution = factor(str_sub(H[[i]]$MutationType, 3, -3), levels = unique(str_sub(H[[i]]$MutationType, 3, -3)))
  H[[i]] = pivot_longer(H[[i]], colnames(H[[i]])[1]:colnames(H[[i]])[sig])
  colnames(H[[i]])[c(3,4)] = c("Signature", "Intensity")
  H[[i]] = H[[i]] %>% arrange(Signature,Substitution)
  H[[i]]$MutationType <- factor(H[[i]]$MutationType, levels = unique(H[[i]]$MutationType))
  p[[idx]] <- ggplot(H[[i]]) + geom_col(aes(x = MutationType, y = Intensity, fill = Substitution))  + theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 40),
          axis.title = element_text(size = 30),
          axis.text.y = element_text(size = 20),
          strip.text = element_text(size = 30),  
          legend.text = element_text(size = 30),  
          legend.title = element_text(size = 30),  
          strip.background = element_rect(colour="white", fill="white") ) + ylim(0,0.4) +
    facet_wrap(vars(Signature), nrow = 4, strip.position = "top", 
               scales="free_x") + ggtitle(modelname[idx])
  idx = idx + 1
}

ggarrange(p[[1]],p[[2]], p[[3]], ncol = 3, common.legend = T, legend = "bottom")

