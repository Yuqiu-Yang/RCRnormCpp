// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//' @title
//' Inv-gamma
//' @description
//' Generate inverse gamma 
//' @export
// [[Rcpp::export]]

NumericVector rinvGamma(int n, double shape, double rate)
{
  NumericVector temp = 1/rgamma(n, shape, 1/rate);
  return temp;
}