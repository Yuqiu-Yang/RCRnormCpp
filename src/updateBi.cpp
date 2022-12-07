// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector findB1(NumericMatrix neg_dat, NumericVector ai,
                     double cc, NumericVector d_neg,
                     int n_patient)
{
  NumericVector result(n_patient);
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = cc * sum(neg_dat(_,i) - ai(i) - d_neg);  
  }
  return result;
}

// [[Rcpp::export]]
NumericVector findB2(NumericMatrix pos_dat, NumericVector ai,
                     NumericVector pos_conc, NumericVector d_pos,
                     int n_patient)
{
  NumericVector result(n_patient);
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = sum((pos_dat(_,i) - ai(i)  - d_pos) * pos_conc);  
  }
  return result;
} 


// [[Rcpp::export]]
NumericVector findB3(NumericMatrix hk_dat, NumericVector ai,
                     NumericVector phi, NumericVector kappa_hk,
                     NumericVector d_hk, 
                     int n_patient, int n_hk)
{
  NumericVector result(n_patient);
  NumericVector temp(n_hk);
  for(int i = 0; i < n_patient; i++)
  {
    temp = kappa_hk[Range(i*n_hk, (i+1)*n_hk - 1)];
    result(i) = sum((phi(i) + temp) * (hk_dat(_,i) - ai(i) - d_hk));  
  }
  return result;
}

// [[Rcpp::export]]
NumericVector findB4(NumericMatrix reg_dat, NumericVector ai,
                     NumericVector phi, NumericVector kappa_reg,
                     NumericVector d_reg, 
                     int n_patient, int n_reg)
{
  NumericVector result(n_patient);
  NumericVector temp(n_reg);
  for(int i = 0; i < n_patient; i++)
  {
    temp = kappa_reg[Range(i*n_reg, (i+1)*n_reg - 1)];
    result(i) = sum((phi(i) + temp) * (reg_dat(_,i) - ai(i) - d_reg));  
  }
  return result;
}


// [[Rcpp::export]]
NumericVector updateBi(List dat, 
                       NumericVector ai,
                       double cc, NumericVector d_neg,
                       NumericVector pos_conc, NumericVector d_pos,
                       NumericVector phi, NumericVector kappa_hk,
                       NumericVector d_hk, 
                       NumericVector kappa_reg,
                       NumericVector d_reg,
                       double mu_b,
                       double sigma2e_neg,
                       double sigma2e,
                       double sigma2b)
{
  
  int n_neg = d_neg.length();
  int n_pos = d_pos.length();
  int n_hk = d_hk.length();
  int n_reg = d_reg.length();
  int n_patient = ai.length();
  
  NumericVector result(n_patient);
  NumericVector B1 = findB1(dat["neg_dat"], ai,
                            cc, d_neg,
                            n_patient);
  NumericVector B2 = findB2(dat["pos_dat"], ai,
                            pos_conc, d_pos,
                            n_patient);
  NumericVector B3 = findB3(dat["hk_dat"], ai,
                            phi, kappa_hk,
                            d_hk, 
                            n_patient, n_hk);
  NumericVector B4 = findB4(dat["reg_dat"], ai,
                            phi, kappa_reg,
                            d_reg, 
                            n_patient, n_reg);
  NumericVector num = B1/sigma2e_neg + (B2 + B3 + B4)/sigma2e + mu_b/sigma2b;
  
  NumericVector temp_hk(n_patient);
  NumericVector temp_reg(n_patient);
  NumericVector temp(n_hk);
  NumericVector temp1(n_reg);
  for(int i = 0; i< n_patient; i++)
  {
    temp = kappa_hk[Range(i*n_hk, (i+1)*n_hk - 1)];
    temp1 = kappa_reg[Range(i*n_reg, (i+1)*n_reg - 1)];
    temp_hk(i) = sum(pow((phi(i) + temp),2));  
    temp_reg(i) = sum(pow((phi(i) + temp1),2));  
  }
  
  NumericVector den = (n_neg * pow(cc, 2))/sigma2e_neg + 
    (sum(pow(pos_conc,2)) + temp_hk + temp_reg)/sigma2e + 1/sigma2b;
  
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = rnorm(1, (num(i))/(den(i)), sqrt(1/(den(i))))[0];
  }
  return result;
}




