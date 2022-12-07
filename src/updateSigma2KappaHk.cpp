// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2KappaHk(List dat,
                           NumericVector kappa_hk,
                           NumericVector lambda_hk,
                           List priors)
{
  NumericMatrix hk_dat = dat["hk_dat"];
  int n_hk = hk_dat.nrow();
  int n_patient = hk_dat.ncol();
  double result;
  double s = 0;
  
  for(int i = 0; i < n_patient; i++)
  {
    for(int h = 0; h < n_hk; h++)
    {
      s += pow(kappa_hk(i * n_hk + h) - lambda_hk(h),2);
    }
  }
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + (n_hk * n_patient)/2, 
                     v + (s)/2)[0];
  return result;
}