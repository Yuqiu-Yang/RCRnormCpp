// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector updateDneg(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         double cc,
                         double sigma2e_neg,
                         double sigma2d_neg)
{
  NumericMatrix neg_dat = dat["neg_dat"];
  int n_patient = ai.length();
  int n_neg = neg_dat.nrow();
  
  NumericVector result(n_neg);
  
  double s;
  double num;
  double den = n_patient/sigma2e_neg + 1/sigma2d_neg;
  
  for(int n=0; n<n_neg; n++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += neg_dat(n, i) - ai(i) - cc * bi(i);
    }
    num = s/sigma2e_neg;
    result(n) = rnorm(1, num/den, sqrt(1/den))[0];
  }
  return result;
}