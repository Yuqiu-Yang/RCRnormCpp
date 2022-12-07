// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2a(List dat,
                     NumericVector ai,
                     double mu_a,
                     List priors)
{
  int n_patient = ai.length();
  double result;
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + n_patient/2, 
                     v + sum(pow(ai-mu_a,2))/2)[0];
  return result;
}