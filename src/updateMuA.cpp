// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

double updateMuA(List dat,
                 NumericVector ai,
                 double sigma2a,
                 List priors)
{
  int n_patient = ai.length();
  double sigma2_mu_a = priors["sigma2_mu_a"];
  double mu_a_hat = priors["mu_a_hat"];
  
  double result;
  double den = n_patient/sigma2a + 1 / sigma2_mu_a;
  double num = sum(ai)/sigma2a + mu_a_hat/sigma2_mu_a;
  result = rnorm(1, num/den, sqrt(1/den))[0];
  return result;
}