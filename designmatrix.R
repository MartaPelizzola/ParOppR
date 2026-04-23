##--------------------------------------------------------
## Factors and design matrices
##--------------------------------------------------------
getdesignMat = function(data3 = NULL,data5 = NULL,data7 = NULL){
  
  designMatres = list()
  
  if(!is.null(data3)){
  ## tri-nucleotide 
  l1 = factor(substr(colnames(data3), start = 1, stop = 1))
  m = factor(substr(colnames(data3), start = 3,stop = 5))
  r1 = factor(substr(colnames(data3), start = 7, stop = 7))
  
  
  Mmono3 = model.matrix(~0+m+l1+r1)
  Mdi3 = model.matrix(~0+m*l1+m*r1)
  Mfull3 = model.matrix(~0+m*l1*r1)
  designMatres[["data3"]] = list(Mmono3 = Mmono3,Mfull3 = Mfull3, Mdi3 = Mdi3)
  }
  
  if(!is.null(data5)){
  ## penta-nucleotide 
  l2 = factor(substr(colnames(data5), start = 1, stop = 1))
  l1 = factor(substr(colnames(data5), start = 2, stop = 2))
  m = factor(substr(colnames(data5), start = 4,stop = 6))
  r1 = factor(substr(colnames(data5), start = 8, stop = 8))
  r2 = factor(substr(colnames(data5), start = 9, stop = 9))
  
  Mmono5 = model.matrix(~0+m+l1+r1+l2+r2)
  Mdi5 = model.matrix(~0+m*l1+m*r1+l1*l2+r1*r2)
  Mfull5 = model.matrix(~0+m*l1*r1*l2*r2)
  designMatres[["data5"]] = list(Mmono5 = Mmono5,Mfull5 = Mfull5, Mdi5 = Mdi5)
  }
  
  if(!is.null(data7)){
  ## septem-nucleotide
  l3 = factor(substr(colnames(data7), start = 1, stop = 1))
  l2 = factor(substr(colnames(data7), start = 2, stop = 2))
  l1 = factor(substr(colnames(data7), start = 3, stop = 3))
  m = factor(substr(colnames(data7), start = 5,stop = 7))
  r1 = factor(substr(colnames(data7), start = 9, stop = 9))
  r2 = factor(substr(colnames(data7), start = 10, stop = 10))
  r3 = factor(substr(colnames(data7), start = 11, stop = 11))
  
  Mmono7 = model.matrix(~0+m+l1+r1+l2+r2+l3+r3)
  Mdi7 = model.matrix(~0+m*l1+m*r1+l1*l2+r1*r2 + l2*l3 + r2*r3)
  Mfull7 = model.matrix(~0+m*l1*r1*l2*r2*l3*r3)
  designMatres[["data7"]] = list(Mmono7 = Mmono7,Mfull7 = Mfull7, Mdi7 = Mdi7)
  }
  return(designMatres)
}
