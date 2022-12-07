// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

double updateMuB(List dat,
                 NumericVector bi,
                 double sigma2b,
                 List priors)
{
  int n_patient = bi.length();
  double result;
  
  double sigma2_mu_b = priors["sigma2_mu_b"];
  double mu_b_hat = priors["mu_b_hat"];
  
  double den = n_patient/sigma2b + 1 / sigma2_mu_b;
  double num = sum(bi)/sigma2b + mu_b_hat/sigma2_mu_b;
  result = rnorm(1, num/den, sqrt(1/den))[0];
  return result;
}