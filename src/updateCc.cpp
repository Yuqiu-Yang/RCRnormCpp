// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateCc(List dat,
                NumericVector ai,
                NumericVector bi,
                NumericVector d_neg,
                double sigma2e_neg,
                List priors)
{
  NumericMatrix neg_dat = dat["neg_dat"];
  int n_patient = ai.length();
  int n_neg = neg_dat.nrow();
  double cc_lb = priors["cc_lb"];
  double cc_ub = priors["cc_ub"];
  double result;
  
  double s = 0;
  double den = n_neg * sum(pow(bi,2));
  for(int i=0; i< n_patient; i++)
  {
    for(int n=0; n<n_neg; n++)
    {
      s += bi(i) * (neg_dat(n, i) - ai(i) - d_neg(n));
    }
  }
  
  result = rtruncatedNorm(1, cc_lb, cc_ub,
                          s/den, sqrt(sigma2e_neg/den))[0];
  return result;
}
