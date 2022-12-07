// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector updateLambdaReg(List dat,
                              NumericVector kappa_reg,
                             double sigma2kappa_reg,
                             List priors)
{
  
  NumericMatrix reg_dat = dat["reg_dat"];
  int n_patient = reg_dat.ncol();
  int n_reg = reg_dat.nrow();
  NumericVector lambda_reg_lb = priors["lambda_reg_lb"];
  NumericVector lambda_reg_ub = priors["lambda_reg_ub"];
  
  NumericVector result(n_reg);
  double s;
  for(int r=0; r < n_reg; r++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += kappa_reg(i*n_reg + r);
    }
    result(r) = rtruncatedNorm(1, lambda_reg_lb(r), lambda_reg_ub(r),
           s/n_patient, sqrt(sigma2kappa_reg/n_patient))[0];
    
  }
  
  return result;
}
