#include <RcppArmadillo.h>
#include <cmath>        // std::abs
#include <tuple>
#include <iostream>
// #include <gperftools/profiler.h>
// "gperftools/profiler.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
double error(arma::colvec y, arma::colvec mu) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += mu[i];
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) - y[i] + mu[i];
    }
  }
  return sum;
}

// [[Rcpp::depends(RcppArmadillo)]]
double error_nb(arma::colvec y, arma::colvec mu, arma::colvec alpha) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++){
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += -(y[i] + alpha[i]) * (log(alpha[i] + y[i]) - log(alpha[i] + mu[i]));
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) -(y[i] + alpha[i]) * (log(alpha[i] + y[i]) - log(alpha[i] + mu[i]));
    }
  }
  return sum;
}

// [[Rcpp::export]]
arma::rowvec glmUpdate(arma::colvec y, arma::mat X, int maxIter = 50, double epsilon=1e-8) {
  if (y.size() <= X.n_cols) {
    //Rcout << "Iterations:";
    return y.as_row();
  }
  //arma::colvec noise = arma::randu(y.size());
  arma::colvec mu = y + 0.1;
  double old = error(y, mu);
  for (int i=0; i<maxIter; i++) {
    arma::colvec z = arma::log(mu) + (y-mu) / mu;
    arma::colvec w = arma::sqrt(mu);
    arma::colvec coef = arma::solve(X.each_col() % w, z % w);
    arma::colvec eta = X * coef;
    mu = arma::exp(eta);
    double newDev = error(y, mu);
    if (2 * std::abs(old - newDev)/(0.1 + std::abs(2*newDev)) < epsilon) {
      break;
    }
    old = newDev;
  }
  return mu.as_row();
}

std::tuple<arma::mat, arma::mat, double> nmf1glm(arma::mat data, arma::mat param[], int rank, arma::colvec alpha, int iter = 5000, arma::rowvec opp = Rcpp::NumericVector::create()) {
  
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  if (opp.size() == 0){
    opp = arma::ones(mutTypes).t();
  }
  
  arma::mat exposures(genomes, rank, arma::fill::randu);
  arma::mat signatures(rank, mutTypes, arma::fill::randu);
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha;
  
  arma::mat estimate = exposures * signatures;
  estimate = estimate.each_row() % opp;
  
  arma::mat sigopp = signatures.each_row() % opp;
  
  double error;
  
  for(int t = 0; t < iter; t++){
    exposures = exposures % (((data/estimate) * arma::trans(sigopp))/(((alphamat + data)/(alphamat + estimate))* arma::trans(sigopp)));
    
    estimate = exposures * signatures;
    estimate = estimate.each_row() % opp;

    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    for(int row = 0; row<rank; row++) {
       signatures.row(row) = glmUpdate(arma::trans(signatures.row(row)), param[row]);
    }
    
    exposures.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    signatures.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    estimate = exposures * signatures;
    estimate = estimate.each_row() % opp;
    sigopp = signatures.each_row() % opp;
    
    if(t - floor(t/10)*10 == 0 && t < 50){
      Rcout << error_nb(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
      Rcout << "--";
    }
    
  }
  
  error = error_nb(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
  //Rcout << "end value:"<<error << "\n";
  
  return {exposures, signatures, error};
}

// [[Rcpp::export]]
List nmfprmnb(arma::mat data, List param, int rank, arma::colvec alpha, arma::rowvec opp = Rcpp::NumericVector::create(), int maxiter = 10000, double tol = 1e-8, int initial = 100, int smallIter = 500) {
  
  int genomes = data.n_rows;
  int mutTypes = data.n_cols;
  
  if (opp.size() == 0){
    opp = arma::ones(mutTypes).t();
  }
  
  arma::mat paramArray[rank];
  for(int i=0; i<rank; i++) {
    NumericMatrix x = param[i];
    paramArray[i] = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  }  
  
  auto res = nmf1glm(data, paramArray, rank, alpha, smallIter, opp);
  auto exposures = std::get<0>(res);
  auto signatures = std::get<1>(res);
  auto errorValue = std::get<2>(res);
  
  for(int i = 1; i < initial; i++){
    auto res = nmf1glm(data, paramArray, rank, alpha, smallIter, opp);
    auto errorNew = std::get<2>(res);
    
    if(errorNew < errorValue){
      errorValue = errorNew;
      exposures = std::get<0>(res);
      signatures = std::get<1>(res);
    }
  }
  
  arma::mat alphamat(genomes, mutTypes, arma::fill::zeros);
  alphamat.each_col() = alpha; 
  arma::mat estimate = exposures * signatures;
  arma::mat sigopp = signatures.each_row() % opp;
  
  double errorOld = error_nb(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
  double errorgkl = error(arma::vectorise(data),arma::vectorise(estimate));
  double errorNew = 2*errorOld;
  double it = 0;
  arma::vec nbobjvalues(maxiter);
  
  for(int t = 0; t < maxiter; t++){
    exposures = exposures % (((data/estimate) * arma::trans(sigopp))/(((alphamat + data)/(alphamat + estimate))* arma::trans(sigopp)));
    
    estimate = exposures * signatures;
    estimate = estimate.each_row() % opp;
    
    signatures = signatures % ((arma::trans(exposures) * (data/estimate))/(arma::trans(exposures)*((alphamat + data)/(alphamat + estimate))));
    
    for(int row = 0; row<rank; row++) {
      signatures.row(row) = glmUpdate(arma::trans(signatures.row(row)), param[row]);
    }
    
    exposures.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    signatures.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    estimate = exposures * signatures;
    estimate = estimate.each_row() % opp;
    sigopp = signatures.each_row() % opp;
    
    if(t - floor(t/10)*10 == 0){
      errorNew = error_nb(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
      
      if (2*std::abs(errorOld - errorNew)/(0.1 + std::abs(2*errorNew)) < tol){
        Rcout << "Iterations:";
        Rcout << t;
        Rcout << "\n";
        Rcout << "end value:"<< errorNew << "\n";
        errorgkl = error(arma::vectorise(data),arma::vectorise(estimate));
        break;
      }
      errorOld = errorNew;
    }
    it = t;
    nbobjvalues.at(t) = errorOld;
  }
  
  arma::colvec rsum = sum(signatures,1);
  exposures = exposures.each_row() % arma::trans(rsum);
  signatures = signatures.each_col() / rsum;
  
  errorgkl = error(arma::vectorise(data),arma::vectorise(estimate));
  
  List output = List::create(Named("exposures") = exposures,
                             Named("signatures") = signatures,
                             Named("gkl") = errorgkl,
                             Named("error") = errorNew,
                             Named("iter") = it,
                             Named("alpha") = alpha,
                             Named("nbobjvalues") = nbobjvalues);
  return output;
}
