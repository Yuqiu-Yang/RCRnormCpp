// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector updateKappaHk(List dat,
                            NumericVector ai,
                            NumericVector bi,
                            NumericVector phi,
                            NumericVector d_hk, 
                            NumericVector lambda_hk,
                            double sigma2e,
                            double sigma2kappa_hk)
{
  int n_hk = d_hk.length();
  int n_patient = ai.length();
  
  NumericVector result(n_patient * n_hk);
  
  NumericMatrix hk_dat = dat["hk_dat"];
  
  int counter = 0;
  
  double den = 0;
  double num = 0;
  for(int i = 0; i < n_patient; i++)
  {
    den = pow(bi(i),2)/sigma2e + 1/sigma2kappa_hk;
    for(int h = 0; h < n_hk; h++)
    {
      num = (bi(i) * (hk_dat(h,i) - ai(i) - bi(i) * phi(i) - d_hk(h)))/sigma2e + lambda_hk(h)/sigma2kappa_hk;
      result(counter) = rnorm(1, 
             num/den,
             sqrt(1/den))[0];
      counter++;
    }
  }
  
  return result;
}


