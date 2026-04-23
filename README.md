# ParOppR
Code to estimate mutational signatures from mutational counts data of cancer patients with the following extensions: (a) incorporating opportunities to the analysis, (b) allowing for extended sequence contexts, (c) using the Negative Binomial model, and (d) parametrizing the signatures.

The data and codes here have been used to create the results seen in the paper 'Integrating opportunities and parametrized signatures for improved mutational processes estimation in extended sequence contexts' by Ragnhild Laursen\*, Marta Pelizzola\*, Lasse Marretty and Asger Hobolth.

The code files ParOppR_mainFunction.R and NegBinNMF_par_opp.cpp are the essential code used to run negative binomial non-negative matrix factorization with parametrization of the signatures and inclusion of mutational opportunities with various different mutational context sizes. 

## How to apply negative binomial non-negative matrix factorization with parametrization of the signatures and inclusion of mutational opportunities
The file ParOppR_mainFunction.R contains the full commented script we used to apply all models presented in the manuscript with information about how to select a given parametrization of the signatures, how is the mutational context size accounted for in the algorithm, and how to include the opportunities in the model.  
Here, you can find a small tutorial where a data set of 21 breast cancer patients is used with simulated opportunities:

1. Download NegBinNMF_par_opp.cpp and designmatrix.R from this repository. 
2. Load both in R:

```r
# libraries used
library(Rcpp)
library(RcppArmadillo)
library(SigMoS) #This contains all support function to apply negative binomial NMF and can be found at https://github.com/MartaPelizzola/SigMoS

# loading the function
sourceCpp("NegBinNMF_par_opp.cpp")
source("designmatrix.R")
```

3. Load your own data set or download the example data from this repository and load it:
```r
data <- load("BRCA21.RData")
```

4. The main function to run negative binomial parametrized NMF with opportunities is 'nmfprmnb' from the cpp file. This takes the dollowing input variables:
  - 'data' = a matrix of count data where the rows correspond to the patients and the columns to the mutation types;
  - 'param' = the design matrix describing the chosen parametrization. The code 'designmatrix.R' in this repository creates the design matrices used in the manuscript from any given data set. For further details on how to define a design matrix see: https://github.com/ragnhildlaursen/paramNMF_ms;
  - 'rank' = the number of signatures to be estimated;
  - 'alpha' = the dispersion parameter for negative binomial NMF. For further details see https://github.com/MartaPelizzola/SigMoS;
  - 'opp' = a numeric vector of opportunities of length equal to ncol(data). Can be empty if running the model without opportunities;
  - 'maxiter' = maximum number of iterations for convergence;
  - 'tol' = tolerance value used for checking convergence;
  - 'inital' = number of initializations;
  - 'smallIter' = number of short-run iterations to test various initializations.

5. In order to run it then we need a vector of opportunities (if desired), the design matrix and alpha:
```r
opp_vec <- rexp(96) # sample mutational opportunities from an exponentail distribution
alpha <- alphaNR(data, 3, TRUE) #estimate alpha using SigMoS
dm <- getdesignMat(data3 = data) #use the function from designmatrix.R to get generate design matrices for the tri-nucleotide context ('data3')
 
```

6. Now, as an example, to apply negative binomial NMF with the interaction parametrization with opportunities and 3 signatures, the following code needs to be used:
```r
design_mat_i <- dm$data3$Mdi3 # extract the design matrix corresponding to the interaction model
out <- nmfprmnb(data, param = rep(list(design_mat_i),3), rank = 3, alpha = alpha,
                 opp = opp_vec, initial = 5, smallIter = 100, maxiter = 10000) #run NMF
```
Note that this is on a tri-nucleotide context as the example data contain mutational counts defined on this context and thus have 96 mutation types.

7. Similarly, to run the same model with the additive parametrization we need:
```r
design_mat_2_a <- dm$data3$Mmono3 # extract the design matrix corresponding to the additive model
out <- nmfprmnb(data, param = rep(list(design_mat_a),3), rank = 3, alpha = alpha,
                 opp = opp_vec, initial = 5, smallIter = 100, maxiter = 10000) #run NMF
```
8. And lastly, to run the same model with the standard parametrization we need:
```r
design_mat_s <- dm$data3$Mfull3 # extract the design matrix corresponding to the standard model
out <- nmfprmnb(data, param = rep(list(design_mat_s),3), rank = 3, alpha = alpha,
                 opp = opp_vec, initial = 5, smallIter = 100, maxiter = 10000) #run NMF
```

9. All other models included in this functions are illustrated in ParOppR_mainFunction.R

## Code for plotting
The code we used to create the figures in our manuscript can be found in the folder 'Figures'.

## Data availability 
The breast and liver cancer data used for the results in this manuscript are available at: https://www.synapse.org/#!Synapse:syn11726620.

## Questions?
If you have any question regarding how to use this code you can contact:
- marta@math.au.dk
- ragnhild@clin.au.dk