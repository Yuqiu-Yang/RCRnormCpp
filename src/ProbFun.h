# ifndef PROBFUN_H
# define PROBFUN_H
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
Rcpp::NumericVector rinvGamma(int n, double shape, double rate);
NumericVector rtruncatedNorm(int n, double lb, double ub, 
                             double mu, double sd);


#endif