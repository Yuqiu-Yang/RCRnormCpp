// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector updateKappaReg(List dat,
                            NumericVector ai,
                            NumericVector bi,
                            NumericVector phi,
                            NumericVector d_reg, 
                            NumericVector lambda_reg,
                            double sigma2e,
                            double sigma2kappa_reg)
{
  int n_reg = d_reg.length();
  int n_patient = ai.length();
  NumericVector result(n_patient * n_reg);
  
  NumericMatrix reg_dat = dat["reg_dat"];
  
  int counter = 0;
  
  double den = 0;
  double num = 0;
  for(int i = 0; i < n_patient; i++)
  {
    den = pow(bi(i),2)/sigma2e + 1/sigma2kappa_reg;
    for(int r = 0; r < n_reg; r++)
    {
      num = (bi(i) * (reg_dat(r,i) - ai(i) - bi(i) * phi(i) - d_reg(r)))/sigma2e + lambda_reg(r)/sigma2kappa_reg;
      result(counter) = rnorm(1, 
             num/den,
             sqrt(1/den))[0];
      counter++;
    }
  }
  
  return result;
}


