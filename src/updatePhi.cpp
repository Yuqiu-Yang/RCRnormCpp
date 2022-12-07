// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector updatePhi(List dat, 
                        NumericVector ai,
                        NumericVector bi,
                        NumericVector kappa_hk,
                        NumericVector d_hk, 
                        NumericVector kappa_reg,
                        NumericVector d_reg,
                        double sigma2e,
                        List priors)
{
  int n_hk = d_hk.length();
  int n_reg = d_reg.length();
  int n_patient = ai.length();
  
  
  NumericVector result(n_patient);
  
  NumericMatrix hk_dat = dat["hk_dat"];
  NumericMatrix reg_dat = dat["reg_dat"];
  
  NumericVector temp_hk(n_patient);
  NumericVector temp_reg(n_patient);
  NumericVector temp(n_hk);
  NumericVector temp1(n_reg);
  
  NumericVector phi_lb = priors["phi_lb"];
  NumericVector phi_ub = priors["phi_ub"];
  for(int i = 0; i< n_patient; i++)
  {
    temp = kappa_hk[Range(i*n_hk, (i+1)*n_hk - 1)];
    temp1 = kappa_reg[Range(i*n_reg, (i+1)*n_reg - 1)];
    temp_hk(i) = sum(hk_dat(_,i) - ai(i) - bi(i) * temp - d_hk);  
    temp_reg(i) = sum(reg_dat(_,i) - ai(i) - bi(i) * temp1 - d_reg);
  }
  
  for(int i = 0; i < n_patient - 1; i++)
  {
    result(i) = rtruncatedNorm(1, phi_lb(i), phi_ub(i),
           (temp_hk(i) + temp_reg(i))/((n_hk + n_reg) * bi(i)),
           sqrt(sigma2e/((n_hk + n_reg) * pow(bi(i),2))))[0];
  }
  result(n_patient-1) = -sum(result[Rcpp::Range(0, n_patient-2)]);
  return result;
}