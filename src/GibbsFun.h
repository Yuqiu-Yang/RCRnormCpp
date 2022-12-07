# ifndef GIBBSFUN_H
# define GIBBSFUN_H
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
List linearRegression(const arma::mat& X, const arma::colvec& y);

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
                       double sigma2a);
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
                       double sigma2b);

NumericVector updatePhi(List dat, 
                        NumericVector ai,
                        NumericVector bi,
                        NumericVector kappa_hk,
                        NumericVector d_hk, 
                        NumericVector kappa_reg,
                        NumericVector d_reg,
                        double sigma2e,
                        List priors);
NumericVector updateKappaHk(List dat,
                            NumericVector ai,
                            NumericVector bi,
                            NumericVector phi,
                            NumericVector d_hk, 
                            NumericVector lambda_hk,
                            double sigma2e,
                            double sigma2kappa_hk);
NumericVector updateKappaReg(List dat,
                             NumericVector ai,
                             NumericVector bi,
                             NumericVector phi,
                             NumericVector d_reg, 
                             NumericVector lambda_reg,
                             double sigma2e,
                             double sigma2kappa_reg);
NumericVector updateLambdaHk(List dat,
                             NumericVector kappa_hk,
                             double sigma2kappa_hk,
                             List priors);
NumericVector updateLambdaReg(List dat,
                              NumericVector kappa_reg,
                              double sigma2kappa_reg,
                              List priors);
NumericVector updateDneg(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         double cc,
                         double sigma2e_neg,
                         double sigma2d_neg);
NumericVector updateDpos(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         NumericVector pos_conc,
                         double sigma2e,
                         double sigma2d);
NumericVector updateDhk(List dat, 
                        NumericVector ai,
                        NumericVector bi,
                        NumericVector phi,
                        NumericVector kappa_hk,
                        double sigma2e,
                        double sigma2d);
NumericVector updateDreg(List dat, 
                         NumericVector ai,
                         NumericVector bi,
                         NumericVector phi,
                         NumericVector kappa_reg,
                         double sigma2e,
                         double sigma2d);
double updateCc(List dat,
                NumericVector ai,
                NumericVector bi,
                NumericVector d_neg,
                double sigma2e_neg,
                List priors);

double updateMuA(List dat,
                 NumericVector ai,
                 double sigma2a,
                 List priors);

double updateMuB(List dat,
                 NumericVector bi,
                 double sigma2b,
                 List priors);

double updateSigma2a(List dat,
                     NumericVector ai,
                     double mu_a,
                     List priors);

double updateSigma2b(List dat,
                     NumericVector bi,
                     double mu_b,
                     List priors);

double updateSigma2Dneg(List dat,
                        NumericVector d_neg,
                        List priors);
double updateSigma2D(List dat,
                     NumericVector d_pos,
                     NumericVector d_hk,
                     NumericVector d_reg,
                     List priors);

double updateSigma2Eneg(List dat,
                        NumericVector ai,
                        NumericVector bi,
                        double cc,
                        NumericVector d_neg,
                        List priors);

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
                     List priors);

double updateSigma2KappaHk(List dat,
                           NumericVector kappa_hk,
                           NumericVector lambda_hk,
                           List priors);

double updateSigma2KappaReg(List dat,
                            NumericVector kappa_reg,
                            NumericVector lambda_reg,
                            List priors);

#endif
