// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector updateDhk(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         NumericVector phi,
                         NumericVector kappa_hk,
                         double sigma2e,
                         double sigma2d)
{
  NumericMatrix hk_dat = dat["hk_dat"];
  int n_patient = ai.length();
  int n_hk = hk_dat.nrow();
  NumericVector result(n_hk);
  
  
  double s;
  double num;
  double den = n_patient/sigma2e + 1/sigma2d;
  
  for(int h=0; h<n_hk; h++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += hk_dat(h, i) - ai(i) - bi(i) * (phi(i) + kappa_hk(i * n_hk + h));
    }
    num = s/sigma2e;
    result(h) = rnorm(1, num/den, sqrt(1/den))[0];
  }
  return result;
}