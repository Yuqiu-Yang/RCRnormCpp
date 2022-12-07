// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ProbFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
double updateSigma2Eneg(List dat,
                        NumericVector ai,
                        NumericVector bi,
                        double cc,
                        NumericVector d_neg,
                        List priors)
{
  NumericMatrix neg_dat = dat["neg_dat"];
  int n_neg = d_neg.length();
  int n_patient = ai.length();
  double result;
  
  double s = 0;
  
  for(int i = 0; i < n_patient; i++)
  {
    for(int n = 0; n < n_neg; n++)
    {
      s += pow(neg_dat(n, i) - ai(i) - bi(i) * cc - d_neg(n),2);
    }
  }
  double u = priors["u"];
  double v = priors["v"];
  result = rinvGamma(1,u + (n_neg * n_patient)/2, 
                     v + (s)/2)[0];
  return result;
}