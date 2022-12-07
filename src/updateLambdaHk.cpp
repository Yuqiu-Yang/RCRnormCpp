// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector updateLambdaHk(List dat,
                             NumericVector kappa_hk,
                             double sigma2kappa_hk,
                             List priors)
{
  NumericMatrix hk_dat = dat["hk_dat"];
  int n_patient = hk_dat.ncol();
  int n_hk = hk_dat.nrow();
  NumericVector lambda_hk_lb =  priors["lambda_hk_lb"];
  NumericVector lambda_hk_ub =  priors["lambda_hk_ub"];
  
  NumericVector result(n_hk);
  double s;
  for(int h=0; h < n_hk; h++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += kappa_hk(i*n_hk + h);
    }
    
    result(h) = rtruncatedNorm(1, lambda_hk_lb(h), lambda_hk_ub(h),
           s/n_patient, sqrt(sigma2kappa_hk/n_patient))[0];
    
  }
  
  return result;
}
