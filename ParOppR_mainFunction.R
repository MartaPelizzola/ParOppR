rm(list = ls())
# This script contains the code to call the functions to estimate mutational 
# signatures with a negative binomial NMF model with parametrized signatures and
# opportunities. 
# It shows how to apply the main functions of this repository to use the 
# various parametrization (additive, interaction and standard) we propose and 
# how to choose whether to include opportunities in the model.

##------------------------------------
## Setup
##------------------------------------
library(SQUAREM)
library(doParallel)
library(foreach)
library(parallel)
library(stats)
library(Rcpp)
library(RcppArmadillo)
library(SigMoS) #This is loading the functions for Negative Binomial NMF and can be found at https://github.com/MartaPelizzola/SigMoS
sourceCpp("NMFoppNB3.cpp")

# data 
load("livercounts3.RData") # A matrix of mutational counts with nrow = number of patients and ncol= 96
load("livercounts5.RData") # A matrix of mutational counts with nrow = number of patients and ncol= 1536
load("livercounts7.RData") # A matrix of mutational counts with nrow = number of patients and ncol= 24576

load("opportunities.RData") # An R object containing 3 vectors of mutational opportunities (opp3, opp5, opp7) 
                            # of lengths 96, 1536, and 24576 respectively. 

source("designmatrix.R")

designMat = getdesignMat(data3 = livercounts3,data5 = livercounts5,data7 = livercounts7)

##--------------------------------------------------------
## Design matrices
##--------------------------------------------------------

# Additive parametrization for the tri, penta and hepta nucleotide contexts #
Mmono3 = designMat$data3$Mmono3
Mmono5 = designMat$data5$Mmono5
Mmono7 = designMat$data7$Mmono7
# Interaction parametrization for the tri, penta and hepta nucleotide contexts #
Mdi3 = designMat$data3$Mdi3
Mdi5 = designMat$data5$Mdi5
Mdi7 = designMat$data7$Mdi7
# Standard parametrization for the tri, penta and hepta nucleotide contexts #
Mfull3 = designMat$data3$Mfull3
Mfull5 = designMat$data5$Mfull5
Mfull7 = designMat$data7$Mfull7

##-----------------------------------------------------------------
## Estimation of dispersion parameter for Negative Binomial NMF
##-----------------------------------------------------------------
sigvec = c(2,4,6,8,10,12,14,16)
alpha3 <- list()
alpha5 <- list()
alpha7 <- list()

for (k in 1:length(sigvec)){
  rank <- sigvec[k]
  alpha3[[k]] <- alphaNR(livercounts3, k = rank, patient_specific = T)
  alpha5[[k]] <- alphaNR(livercounts5, k = rank, patient_specific = T)
  alpha7[[k]] <- alphaNR(livercounts7, k = rank, patient_specific = T)
}
saveRDS(alpha3, "liver_alpha3.rds")
saveRDS(alpha5,	"liver_alpha5.rds")
saveRDS(alpha7,	"liver_alpha7.rds")

##--------------------------------------------------------
## Run NMF 
##--------------------------------------------------------


init = 3 # Number of initializations
siter = 500 # Number of iterations for all 'init' initializations
miter = 10000 # Maximum number of iterations for the best initializations (if convergence is not reached)
liverres = list()
liverres$noSig = sigvec
idx = 1
 
for(noSig in sigvec){
  # Additive parametrization ('mono') for the tri-nucleotide context (3) without opportunities 
  
  out = nmfprmnb(livercounts3, param = rep(list(Mmono3),noSig), rank = noSig, alpha = alpha3[[idx]],
                 initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono3"]][[idx]] = out
  
  # Additive parametrization ('mono') for the penta-nucleotide context (5) without opportunities
  
  out = nmfprmnb(livercounts5, param = rep(list(Mmono5),noSig), rank = noSig, alpha = alpha5[[idx]],
                 initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono5"]][[idx]] = out
  
  # Additive parametrization ('mono') for the hepta-nucleotide context (7) without opportunities
  
  out = nmfprmnb(livercounts7, param = rep(list(Mmono7),noSig), rank = noSig, alpha = alpha7[[idx]],
                 initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono7"]][[idx]] = out
  
  # Additive parametrization ('mono') for the tri-nucleotide context (3) with opportunities ('opp') 
 
  out = nmfprmnb(livercounts3, param = rep(list(Mmono3),noSig), rank = noSig, alpha = alpha3[[idx]],
                 opp = opp3, initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono3opp"]][[idx]] = out
  
  # Additive parametrization ('mono') for the penta-nucleotide context (5) with opportunities ('opp') 
  
  out = nmfprmnb(livercounts5, param = rep(list(Mmono5),noSig), rank = noSig, alpha = alpha5[[idx]],
               opp = opp5, initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono5opp"]][[idx]] = out
  
  # Additive parametrization ('mono') for the hepta-nucleotide context (7) with opportunities ('opp') 
  
  out = nmfprmnb(livercounts7, param = rep(list(Mmono7),noSig), rank = noSig, alpha = alpha7[[idx]],
                 opp = opp7, initial = init, smallIter = siter, maxiter = miter)
  liverres[["mono7opp"]][[idx]] = out

  
  # Interaction parametrization ('di') for the tri-nucleotide context (3) without opportunities
  
  out = nmfprmnb(livercounts3, param = rep(list(Mdi3),noSig), rank = noSig, alpha = alpha3[[idx]],
               initial = init, smallIter = siter, maxiter = miter)
  liverres[["di3"]][[idx]] = out
  
  # Interaction parametrization ('di') for the penta-nucleotide context (5) without opportunities 
  
  out = nmfprmnb(livercounts5, param = rep(list(Mdi5),noSig), rank = noSig, alpha = alpha5[[idx]],
               initial = init, smallIter = siter, maxiter = miter)
  liverres[["di5"]][[idx]] = out
  
  # Interaction parametrization ('di') for the hepta-nucleotide context (7) without opportunities
  
  out = nmfprmnb(livercounts7, param = rep(list(Mdi7),noSig), rank = noSig, alpha = alpha7[[idx]],
                 initial = init, smallIter = siter, maxiter = miter)
  liverres[["di7"]][[idx]] = out
  
  
  # Interaction parametrization ('di') for the tri-nucleotide context (3) with opportunities ('opp') 
  
  out = nmfprmnb(livercounts3, param = rep(list(Mdi3),noSig), rank = noSig, alpha = alpha3[[idx]],
                 opp = opp3, initial = init, smallIter = siter, maxiter = miter)
  liverres[["di3opp"]][[idx]] = out
  
  # Interaction parametrization ('di') for the penta-nucleotide context (5) with opportunities ('opp') 
  
  out = nmfprmnb(livercounts5, param = rep(list(Mdi5),noSig), rank = noSig, alpha = alpha5[[idx]],
                 opp = opp5, initial = init, smallIter = siter, maxiter = miter)
  liverres[["di5opp"]][[idx]] = out
  
  # Interaction parametrization ('di') for the hepta-nucleotide context (7) with opportunities ('opp') 
  
  out = nmfprmnb(livercounts7, param = rep(list(Mdi7),noSig), rank = noSig, alpha = alpha7[[idx]],
                 opp = opp7, initial = init, smallIter = siter, maxiter = miter)
  liverres[["di7opp"]][[idx]] = out
  
  
  # Standard parametrization ('full') for the tri-nucleotide context (3) 
  
  out = nmfprmnb(livercounts3, param = rep(list(Mfull3),noSig), rank = noSig, alpha = alpha3[[idx]],
                initial = init, smallIter = siter, maxiter = miter)
  liverres[["full3"]][[idx]] = out
  
  # Standard parametrization ('full') for the penta-nucleotide context (5) 
  
  out = nmfprmnb(livercounts5, param = rep(list(Mfull5),noSig), rank = noSig, alpha = alpha5[[idx]],
               initial = init, smallIter = siter, maxiter = miter)
  liverres[["full5"]][[idx]] = out
  
  # Standard parametrization ('full') for the hepta-nucleotide context (7) 
  
  out = nmfprmnb(livercounts7, param = rep(list(Mfull7),noSig), rank = noSig, alpha = alpha7[[idx]],
                 initial = init, smallIter = siter, maxiter = miter)
  liverres[["full7"]][[idx]] = out
   
  
  cat("Finished running", noSig, "signatures \n")
  idx = idx + 1
  
  save(liverres, file = "LIVERres3_5ThreeInitNBwithOPP.RData")
}
