// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2D(List dat,
                     NumericVector d_pos,
                     NumericVector d_hk,
                     NumericVector d_reg,
                     List priors)
{
  int n_pos = d_pos.length();
  int n_hk = d_hk.length();
  int n_reg = d_reg.length();
  double result;
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + (n_pos + n_hk + n_reg)/2, 
                     v + (sum(pow(d_pos,2)) + sum(pow(d_hk,2)) + sum(pow(d_reg,2)))/2)[0];
  return result;
}