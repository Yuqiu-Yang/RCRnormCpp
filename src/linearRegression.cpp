// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//' @title
//' Linear Regression
//' @description
//' Performs basic linear regression
//' @param X explanatory variables (should be a matrix)
//' @param y the response variable 
//' @return
//' A list containing coefficients and residuals
//' @export
// [[Rcpp::export]]

Rcpp::List linearRegression(const arma::mat& X, const arma::colvec& y)
{
  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X * coef;
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("residuals") = resid);
}