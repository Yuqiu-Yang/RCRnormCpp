// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double lbLambdaUb(double lambda, double lb, double ub)
{
  double result = -log(lambda) + pow(lambda, 2)/2 - lambda * lb;
  return result;
}


// [[Rcpp::export]]
double lambdaLbUb(double lambda, double lb, double ub)
{
  double result = -log(lambda) - pow(lb, 2)/2;
  return result;
}


// [[Rcpp::export]]
double lbUbLambda(double lambda, double lb, double ub)
{
  double result = -log(lambda) - pow(ub, 2)/2 + lambda* (ub - lb);
  return result;
}

// [[Rcpp::export]]
NumericVector rtruncatedStandardNormPositive(int n, double lb, double ub)
{
  NumericVector result(n);
  double temp1 = (lb + sqrt(pow(lb, 2) + 4))/2;
  double temp2 = 1/(ub - lb);
  NumericVector temp(3);
  double lambda = 0;
  
  double v = 0;
  double u = 0;
  double ratio = 0;
  
  if((ub <= temp1) & (ub <= temp2))
  {
    temp = NumericVector::create(lbLambdaUb(ub, lb, ub),
                                 lambdaLbUb(lb, lb, ub),
                                 lbUbLambda(temp2, lb, ub));
    lambda = NumericVector::create(ub, lb, temp2)[which_min(temp)];
  }else if((ub <= temp1) & (ub > temp2)){
    temp = NumericVector::create(lbLambdaUb(ub, lb, ub),
                                 lambdaLbUb(lb, lb, ub),
                                 lbUbLambda(ub, lb, ub));
    lambda = NumericVector::create(ub, lb, ub)[which_min(temp)];
  }else if((ub > temp1) & (ub <= temp2)){
    temp = NumericVector::create(lbLambdaUb(temp1, lb, ub),
                                 lambdaLbUb(lb, lb, ub),
                                 lbUbLambda(temp2, lb, ub));
    lambda = NumericVector::create(temp1, lb, temp2)[which_min(temp)];
  }else{
    temp = NumericVector::create(lbLambdaUb(temp1, lb, ub),
                                 lambdaLbUb(lb, lb, ub),
                                 lbUbLambda(ub, lb, ub));  
    lambda = NumericVector::create(temp1, lb, ub)[which_min(temp)];
  }

  double multi_const = (1-exp(-(ub-lb) * lambda))/
    (sqrt(2 * M_PI) * (pnorm(NumericVector::create(ub))[0] - 
      pnorm(NumericVector::create(lb))[0]));
  
  double M = min(temp);
  M += log(multi_const);
  bool accept = false;
  for(int i = 0; i < n; i++)
  {
    accept = false;
    while(! accept)
    {
      v = ub;
      while(v >= ub - lb)
      {
        v = rexp(1, lambda)[0];
      }
      v += lb;
      ratio = (dnorm(NumericVector::create(v))[0])/(pnorm(NumericVector::create(ub))[0] - 
        pnorm(NumericVector::create(lb))[0]);
      ratio /= (lambda * exp(-lambda* (v - lb)))/(1 - exp(-lambda* (ub - lb)));
      u = runif(1)[0];
      if(log(u) < log(ratio) - M)
      {
        result[i] = v;
        accept = true;
      }
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericVector rtruncatedStandardNormPositiveExtreme(int n, double lb, double ub)
{
  NumericVector result(n);
  double x = 0;
  double rho = 0;
  double u = 0;
  for(int i =0; i < n; i++)
  {
    x = ub + 1;
    rho = -1;
    while(u > rho | x > ub)
    {
      x = rexp(1, lb)[0] + lb;
      rho = exp(-0.5 * pow((x - lb) , 2));
      u = runif(1)[0];
    }
    result(i) = x;
  }
  return result;
}


// [[Rcpp::export]]
NumericVector rtruncatedStandardNormNegative(int n, double lb, double ub)
{
  NumericVector temp = rtruncatedStandardNormPositive(n, -ub, -lb);
  return(-temp);
}

// [[Rcpp::export]]
NumericVector rtruncatedStandardNormNegativeExtreme(int n, double lb, double ub)
{
  NumericVector temp = rtruncatedStandardNormPositiveExtreme(n, -ub, -lb);
  return(-temp);
}


// [[Rcpp::export]]
NumericVector rtruncatedStandardNorm(int n, double lb, double ub)
{
  NumericVector temp(n);
  
  double p = pnorm(NumericVector::create(ub))[0] - pnorm(NumericVector::create(0))[0];
  double pos_p = max(NumericVector::create(0, p));
  p = pnorm(NumericVector::create(0))[0] - pnorm(NumericVector::create(lb))[0];
  double neg_p = max(NumericVector::create(0, p));
  pos_p /= (pos_p + neg_p);
  neg_p = 1 - pos_p;
  int pos_n = rbinom(1, n, pos_p)[0];
  int neg_n = n - pos_n;
  
  
  if(neg_p == 0)
  {
    if(lb >= 5)
    {
      temp = rtruncatedStandardNormPositiveExtreme(n, lb, ub);
    }else{
      temp = rtruncatedStandardNormPositive(n, lb, ub);
    }
    
  }else if(pos_p == 0){
    if(ub <= -5)
    {
      temp = rtruncatedStandardNormNegativeExtreme(n, lb, ub);
    }else{
      temp = rtruncatedStandardNormNegative(n, lb, ub);
    }
  }else{
    if(neg_n == 0)
    {
      temp[Rcpp::Range(0, pos_n-1)] = rtruncatedStandardNormPositive(pos_n, 0, ub);
    }else if(pos_n == 0){
      temp[Rcpp::Range(pos_n, n - 1)] = rtruncatedStandardNormNegative(neg_n, lb, 0);
    }else{
      temp[Rcpp::Range(0, pos_n-1)] = rtruncatedStandardNormPositive(pos_n, 0, ub);
      temp[Rcpp::Range(pos_n, n - 1)] = rtruncatedStandardNormNegative(neg_n, lb, 0);
    }
    
  }
  return temp;
}


// [[Rcpp::export]]
NumericVector rtruncatedNorm(int n, double lb, double ub, 
                             double mu, double sd)
{
  NumericVector result(n);
  double lb_standard = (lb - mu)/sd;
  double ub_standard = (ub - mu)/sd;
  
  result = rtruncatedStandardNorm(n, lb_standard, ub_standard);
  result = result * sd;
  result = result + mu;
  return result;
}

