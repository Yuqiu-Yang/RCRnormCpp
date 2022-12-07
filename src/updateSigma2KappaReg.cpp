// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2KappaReg(List dat,
                            NumericVector kappa_reg,
                            NumericVector lambda_reg,
                            List priors)
{
  NumericMatrix reg_dat = dat["reg_dat"];
  int n_reg = reg_dat.nrow();
  int n_patient = reg_dat.ncol();
  
  double result;
  double s = 0;
  
  for(int i = 0; i < n_patient; i++)
  {
    for(int r = 0; r < n_reg; r++)
    {
      s += pow(kappa_reg(i * n_reg + r) - lambda_reg(r),2);
    }
  }
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + (n_reg * n_patient)/2, 
                     v + (s)/2)[0];
  return result;
}