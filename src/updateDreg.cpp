// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector updateDreg(List dat, 
                        NumericVector ai,
                        NumericVector bi,
                        NumericVector phi,
                        NumericVector kappa_reg,
                        double sigma2e,
                        double sigma2d)
{
  NumericMatrix reg_dat = dat["reg_dat"];
  int n_patient = ai.length();
  int n_reg = reg_dat.nrow();
  NumericVector result(n_reg);
  
  
  double s;
  double num;
  double den = n_patient/sigma2e + 1/sigma2d;
  
  for(int r=0; r<n_reg; r++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += reg_dat(r, i) - ai(i) - bi(i) * (phi(i) + kappa_reg(i * n_reg + r));
    }
    num = s/sigma2e;
    result(r) = rnorm(1, num/den, sqrt(1/den))[0];
  }
  return result;
}