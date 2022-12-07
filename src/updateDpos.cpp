// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector updateDpos(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         NumericVector pos_conc,
                         double sigma2e,
                         double sigma2d)
{
  NumericMatrix pos_dat = dat["pos_dat"];
  int n_patient = ai.length();
  int n_pos = pos_dat.nrow();
  NumericVector result(n_pos);
  
  
  double s;
  double num;
  double den = n_patient/sigma2e + 1/sigma2d;
  
  for(int p=0; p<n_pos; p++)
  {
    s = 0;
    for(int i = 0; i < n_patient; i++)
    {
      s += pos_dat(p, i) - ai(i) -  pos_conc(p) * bi(i);
    }
    num = s/sigma2e;
    result(p) = rnorm(1, num/den, sqrt(1/den))[0];
  }
  return result;
}