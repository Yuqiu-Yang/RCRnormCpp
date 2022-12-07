// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector findA1(NumericMatrix neg_dat, NumericVector bi,
                 double cc, NumericVector d_neg,
                 int n_patient)
{
  // Each column of neg_dat should be a patient
  // Each row should be a gene
  // For patient i, A1i is found by sum_{n=1}^N(Y_{ni} - b_i *c - d_n)
  NumericVector result(n_patient);
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = sum(neg_dat(_,i) - bi(i) * cc - d_neg);  
  }
  return result;
}


// [[Rcpp::export]]
NumericVector findA2(NumericMatrix pos_dat, NumericVector bi,
                 NumericVector pos_conc, NumericVector d_pos,
                 int n_patient)
{
  NumericVector result(n_patient);
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = sum(pos_dat(_,i) - bi(i) * pos_conc - d_pos);  
  }
  return result;
}

// [[Rcpp::export]]
NumericVector findA3(NumericMatrix hk_dat, NumericVector bi,
                 NumericVector phi, NumericVector kappa_hk,
                 NumericVector d_hk, 
                 int n_patient, int n_hk)
{
  NumericVector result(n_patient);
  NumericVector temp(n_hk);
  for(int i = 0; i < n_patient; i++)
  {
    temp = kappa_hk[Range(i*n_hk, (i+1)*n_hk - 1)];
    result(i) = sum(hk_dat(_,i) - bi(i) * (phi(i) + temp) - d_hk);  
  }
  return result;
}

// [[Rcpp::export]]
NumericVector findA4(NumericMatrix reg_dat, NumericVector bi,
                 NumericVector phi, NumericVector kappa_reg,
                 NumericVector d_reg, 
                 int n_patient, int n_reg)
{
  NumericVector result(n_patient);
  NumericVector temp(n_reg);
  for(int i = 0; i < n_patient; i++)
  {
    temp = kappa_reg[Range(i*n_reg, (i+1)*n_reg - 1)];
    result(i) = sum(reg_dat(_,i) - bi(i) * (phi(i) + temp) - d_reg);  
  }
  return result;
}


// [[Rcpp::export]]
NumericVector updateAi(List dat, 
                       NumericVector bi,
                       double cc, NumericVector d_neg,
                       NumericVector pos_conc, NumericVector d_pos,
                       NumericVector phi, NumericVector kappa_hk,
                       NumericVector d_hk, 
                       NumericVector kappa_reg,
                       NumericVector d_reg,
                       double mu_a,
                       double sigma2e_neg,
                       double sigma2e,
                       double sigma2a)
{
  int n_neg = d_neg.length();
  int n_pos = d_pos.length();
  int n_hk = d_hk.length();
  int n_reg = d_reg.length();
  int n_patient = bi.length();
  
  
  NumericVector result(n_patient);
  NumericVector A1 = findA1(dat["neg_dat"], bi,
                            cc, d_neg,
                            n_patient);
  
  NumericVector A2 = findA2(dat["pos_dat"], bi,
                       pos_conc, d_pos,
                       n_patient);

  NumericVector A3 = findA3(dat["hk_dat"], bi,
                            phi, kappa_hk,
                            d_hk, 
                            n_patient, n_hk);

  NumericVector A4 = findA4(dat["reg_dat"], bi,
                            phi, kappa_reg,
                            d_reg, 
                            n_patient, n_reg);

  NumericVector num = A1/sigma2e_neg + (A2 + A3 + A4)/sigma2e + mu_a/sigma2a;
  double den = n_neg/sigma2e_neg + (n_pos + n_hk + n_reg)/sigma2e + 1/sigma2a;
  
  for(int i = 0; i < n_patient; i++)
  {
    result(i) = rnorm(1, (num(i))/den, sqrt(1/den))[0];
  }
  return result;
}



