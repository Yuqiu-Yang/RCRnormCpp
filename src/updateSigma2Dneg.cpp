// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2Dneg(List dat,
                        NumericVector d_neg,
                        List priors)
{
  int n_neg = d_neg.length();
  double result;
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + n_neg/2, 
                     v + sum(pow(d_neg,2))/2)[0];
  return result;
}