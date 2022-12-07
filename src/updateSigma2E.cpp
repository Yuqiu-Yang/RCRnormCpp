// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2E(List dat,
                     NumericVector ai,
                     NumericVector bi,
                     NumericVector pos_conc,
                     NumericVector phi,
                     NumericVector kappa_hk,
                     NumericVector kappa_reg,
                     NumericVector d_pos,
                     NumericVector d_hk,
                     NumericVector d_reg,
                     List priors)
{
  NumericMatrix pos_dat = dat["pos_dat"];
  NumericMatrix hk_dat = dat["hk_dat"];
  NumericMatrix reg_dat = dat["reg_dat"];
  int n_pos = d_pos.length();
  int n_hk = d_hk.length();
  int n_reg = d_reg.length();
  int n_patient = pos_dat.ncol();
  double result;
  
  double s1 = 0;
  double s2 = 0;
  double s3 = 0;
  
  for(int i = 0; i < n_patient; i++)
  {
    for(int p = 0; p < n_pos; p++)
    {
      s1 += pow(pos_dat(p, i) - ai(i) - bi(i) * pos_conc(p) - d_pos(p),2);
    }
    for(int h = 0; h < n_hk; h++)
    {
      s2 += pow(hk_dat(h, i) - ai(i) - bi(i) * (phi(i) + kappa_hk(i * n_hk + h)) - d_hk(h),2);
    }
    for(int r = 0; r < n_reg; r++)
    {
      s3 += pow(reg_dat(r, i) - ai(i) - bi(i) * (phi(i) + kappa_reg(i * n_reg + r)) - d_reg(r),2);
    }
  }
  
  double s = s1 + s2 + s3;
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + (n_patient * (n_pos + n_hk + n_reg))/2, 
                     v + (s)/2)[0];
  return result;
}