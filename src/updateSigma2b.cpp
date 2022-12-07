// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2b(List dat,
                     NumericVector bi,
                     double mu_b,
                     List priors)
{
  int n_patient = bi.length();
  double result;
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + n_patient/2, 
                     v + sum(pow(bi-mu_b,2))/2)[0];
  return result;
}