// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "GibbsFun.h"
using namespace Rcpp;
//' @title
//' Gibbs sampler
//' @description
//' One MCMC chain
//' @export

// [[Rcpp::export]]
NumericMatrix generateSampleMatrix(NumericVector p, int n_iter)
{
  NumericMatrix result(n_iter, p.length());
  result(0,_) = p;
  return result;
}


// [[Rcpp::export]]

List sampler(List dat, 
             List priors,
             List initial_values,
             NumericVector pos_conc,
             int n_iter)
{
  // Unpack the variables 
  // NumericMatrix pos_dat = dat["pos_dat"];
  // NumericMatrix neg_dat = dat["neg_dat"];
  // NumericMatrix hk_dat = dat["hk_dat"];
  // NumericMatrix reg_dat = dat["reg_dat"];
  
  // NumericMatrix ai_sample = sample_matrices["ai_sample"];
  // NumericMatrix bi_sample = sample_matrices["bi_sample"];
  // NumericMatrix mu_a_sample = sample_matrices["mu_a_sample"];
  // NumericMatrix mu_b_sample = sample_matrices["mu_b_sample"];
  // NumericMatrix sigma2a_sample = sample_matrices["sigma2a_sample"];
  // NumericMatrix sigma2b_sample = sample_matrices["sigma2b_sample"];
  // NumericMatrix cc_sample = sample_matrices["cc_sample"];
  // NumericMatrix phi_sample = sample_matrices["phi_sample"];
  // NumericMatrix kappa_hk_sample = sample_matrices["kappa_hk_sample"];
  // NumericMatrix kappa_reg_sample = sample_matrices["kappa_reg_sample"];
  // NumericMatrix lambda_hk_sample = sample_matrices["lambda_hk_sample"];
  // NumericMatrix lambda_reg_sample = sample_matrices["lambda_reg_sample"];
  // NumericMatrix sigma2kappa_hk_sample = sample_matrices["sigma2kappa_hk_sample"];
  // NumericMatrix sigma2kappa_reg_sample = sample_matrices["sigma2kappa_reg_sample"];
  // NumericMatrix d_pos_sample = sample_matrices["d_pos_sample"];
  // NumericMatrix d_neg_sample = sample_matrices["d_neg_sample"];
  // NumericMatrix d_hk_sample = sample_matrices["d_hk_sample"];
  // NumericMatrix d_reg_sample = sample_matrices["d_reg_sample"];
  // NumericMatrix sigma2d_sample = sample_matrices["sigma2d_sample"];
  // NumericMatrix sigma2d_neg_sample = sample_matrices["sigma2d_neg_sample"];
  // NumericMatrix sigma2e_sample = sample_matrices["sigma2e_sample"];
  // NumericMatrix sigma2e_neg_sample = sample_matrices["sigma2e_neg_sample"];
  // 
  // double u = priors["u"];
  // double v = priors["v"];
  // double mu_a_hat = priors["mu_a_hat"];
  // double mu_b_hat = priors["mu_b_hat"];
  // double sigma2_mu_a = priors["sigma2_mu_a"];
  // double sigma2_mu_b = priors["sigma2_mu_b"];
  // double lambda_hk_lb = priors["lambda_hk_lb"];
  // double lambda_hk_ub = priors["lambda_hk_ub"]; 
  // double lambda_reg_lb = priors["lambda_reg_lb"];
  // double lambda_reg_ub = priors["lambda_reg_ub"]; 
  // double phi_lb = priors["phi_lb"];
  // double phi_ub = priors["phi_ub"];
  // double cc_lb = priors["cc_lb"];
  // double cc_ub = priors["cc_ub"];
  
  
  
  n_iter++;
  
  NumericMatrix ai_sample = generateSampleMatrix(initial_values["ai"], n_iter);
  NumericMatrix bi_sample = generateSampleMatrix(initial_values["bi"], n_iter);
  NumericMatrix mu_a_sample = generateSampleMatrix(initial_values["mu_a"], n_iter);
  NumericMatrix mu_b_sample = generateSampleMatrix(initial_values["mu_b"], n_iter);
  NumericMatrix sigma2a_sample = generateSampleMatrix(initial_values["sigma2a"], n_iter);
  NumericMatrix sigma2b_sample = generateSampleMatrix(initial_values["sigma2b"], n_iter);
  NumericMatrix cc_sample = generateSampleMatrix(initial_values["cc"], n_iter);
  NumericMatrix phi_sample = generateSampleMatrix(initial_values["phi"], n_iter);
  NumericMatrix kappa_hk_sample = generateSampleMatrix(initial_values["kappa_hk"], n_iter);
  NumericMatrix kappa_reg_sample = generateSampleMatrix(initial_values["kappa_reg"], n_iter);
  NumericMatrix lambda_hk_sample = generateSampleMatrix(initial_values["lambda_hk"], n_iter);
  NumericMatrix lambda_reg_sample = generateSampleMatrix(initial_values["lambda_reg"], n_iter);
  NumericMatrix sigma2kappa_hk_sample = generateSampleMatrix(initial_values["sigma2kappa_hk"], n_iter);
  NumericMatrix sigma2kappa_reg_sample = generateSampleMatrix(initial_values["sigma2kappa_reg"], n_iter);
  NumericMatrix d_pos_sample = generateSampleMatrix(initial_values["d_pos"], n_iter);
  NumericMatrix d_neg_sample = generateSampleMatrix(initial_values["d_neg"], n_iter);
  NumericMatrix d_hk_sample = generateSampleMatrix(initial_values["d_hk"], n_iter);
  NumericMatrix d_reg_sample = generateSampleMatrix(initial_values["d_reg"], n_iter);
  NumericMatrix sigma2d_sample = generateSampleMatrix(initial_values["sigma2d"], n_iter);
  NumericMatrix sigma2d_neg_sample = generateSampleMatrix(initial_values["sigma2d_neg"], n_iter);
  NumericMatrix sigma2e_sample = generateSampleMatrix(initial_values["sigma2e"], n_iter);
  NumericMatrix sigma2e_neg_sample = generateSampleMatrix(initial_values["sigma2e_neg"], n_iter);
  
  NumericVector ai = initial_values["ai"];
  NumericVector bi = initial_values["bi"];
  double  mu_a = initial_values["mu_a"];
  double  mu_b = initial_values["mu_b"];
  double  sigma2a = initial_values["sigma2a"];
  double  sigma2b = initial_values["sigma2b"];
  double  cc = initial_values["cc"];
  NumericVector phi = initial_values["phi"];
  NumericVector kappa_hk = initial_values["kappa_hk"];
  NumericVector kappa_reg = initial_values["kappa_reg"];
  NumericVector lambda_hk = initial_values["lambda_hk"];
  NumericVector lambda_reg = initial_values["lambda_reg"];
  double  sigma2kappa_hk = initial_values["sigma2kappa_hk"];
  double  sigma2kappa_reg = initial_values["sigma2kappa_reg"];
  NumericVector d_pos = initial_values["d_pos"];
  NumericVector d_neg = initial_values["d_neg"];
  NumericVector d_hk = initial_values["d_hk"];
  NumericVector d_reg = initial_values["d_reg"];
  double  sigma2d = initial_values["sigma2d"];
  double  sigma2d_neg = initial_values["sigma2d_neg"];
  double  sigma2e = initial_values["sigma2e"];
  double  sigma2e_neg = initial_values["sigma2e_neg"];

  
  for(int iter = 1; iter < n_iter; iter++)
  {
    /////////////////////
    // update ai
    ////////////////////
    ai = updateAi(dat = dat,
                  bi = bi,
                  cc = cc,
                  d_neg = d_neg,
                  pos_conc = pos_conc,
                  d_pos = d_pos,
                  phi = phi, 
                  kappa_hk = kappa_hk,
                  d_hk = d_hk, 
                  kappa_reg = kappa_reg,
                  d_reg = d_reg,
                  mu_a = mu_a,
                  sigma2e_neg = sigma2e_neg,
                  sigma2e = sigma2e,
                  sigma2a = sigma2a);
    ai_sample(iter, _) = ai;
    //Rcout << "ai" << "\n";
    ////////////////////
    // update bi
    ////////////////////
    bi = updateBi(dat = dat, 
                  ai = ai,
                  cc = cc, 
                  d_neg = d_neg,
                  pos_conc = pos_conc, 
                  d_pos = d_pos,
                  phi = phi, 
                  kappa_hk = kappa_hk,
                  d_hk = d_hk, 
                  kappa_reg = kappa_reg,
                  d_reg = d_reg,
                  mu_b = mu_b,
                  sigma2e_neg = sigma2e_neg,
                  sigma2e = sigma2e,
                  sigma2b = sigma2b);
    bi_sample(iter,_) = bi;
    //Rcout << "bi" << "\n";
    /////////////////////
    // update phi
    ////////////////////
    phi = updatePhi(dat = dat, 
                    ai = ai,
                    bi = bi,
                    kappa_hk = kappa_hk,
                    d_hk = d_hk, 
                    kappa_reg = kappa_reg,
                    d_reg = d_reg,
                    sigma2e = sigma2e,
                    priors = priors);
    phi_sample(iter,_) = phi;
    //Rcout << "phi" << "\n";
    /////////////////////
    // update kappa_hk
    ////////////////////
    kappa_hk = updateKappaHk(dat = dat,
                             ai = ai,
                             bi = bi,
                             phi = phi,
                             d_hk = d_hk, 
                             lambda_hk = lambda_hk,
                             sigma2e = sigma2e,
                             sigma2kappa_hk = sigma2kappa_hk);
    kappa_hk_sample(iter, _) = kappa_hk;
    //Rcout << "khk" << "\n";
    /////////////////////
    // update kappa_reg
    ////////////////////
    kappa_reg = updateKappaReg(dat = dat,
                             ai = ai,
                             bi = bi,
                             phi = phi,
                             d_reg = d_reg, 
                             lambda_reg = lambda_reg,
                             sigma2e = sigma2e,
                             sigma2kappa_reg = sigma2kappa_reg);
    kappa_reg_sample(iter, _) = kappa_reg;
    //Rcout << "kreg" << "\n";
    /////////////////////
    // update lambda_hk
    ////////////////////
    lambda_hk = updateLambdaHk(dat = dat,
                               kappa_hk = kappa_hk,
                               sigma2kappa_hk = sigma2kappa_hk,
                               priors = priors);
    lambda_hk_sample(iter, _) = lambda_hk;
    //Rcout << "lhk" << "\n";
    /////////////////////
    // update lambda_reg
    ////////////////////
    lambda_reg = updateLambdaReg(dat = dat,
                               kappa_reg = kappa_reg,
                               sigma2kappa_reg = sigma2kappa_reg,
                               priors = priors);
    lambda_reg_sample(iter, _) = lambda_reg;
    //Rcout << "lreg" << "\n";
    /////////////////////
    // update d_neg
    ////////////////////
    d_neg = updateDneg(dat = dat,
                       ai = ai,
                       bi = bi,
                       cc = cc,
                       sigma2e_neg = sigma2e_neg,
                       sigma2d_neg = sigma2d_neg);
    d_neg_sample(iter, _) = d_neg;
    //Rcout << "dneg" << "\n";
    /////////////////////
    // update d_pos
    ////////////////////
    d_pos = updateDpos(dat = dat,
                       ai = ai,
                       bi = bi,
                       pos_conc = pos_conc,
                       sigma2e = sigma2e,
                       sigma2d = sigma2d);
    d_pos_sample(iter,_) = d_pos;
    //Rcout << "dpos" << "\n";
    /////////////////////
    // update d_hk
    ////////////////////
    d_hk = updateDhk(dat = dat,
                     ai = ai,
                     bi = bi,
                     phi = phi,
                     kappa_hk = kappa_hk,
                     sigma2e = sigma2e,
                     sigma2d = sigma2d);
    d_hk_sample(iter,_) = d_hk;
    //Rcout << "dhk" << "\n";

    /////////////////////
    // update d_reg
    ////////////////////
    d_reg = updateDreg(dat = dat,
                     ai = ai,
                     bi = bi,
                     phi = phi,
                     kappa_reg = kappa_reg,
                     sigma2e = sigma2e,
                     sigma2d = sigma2d);
    d_reg_sample(iter,_) = d_reg;
    //Rcout << "dreg" << "\n";
    /////////////////////
    // update cc
    ////////////////////

    cc = updateCc(dat = dat,
                  ai = ai,
                  bi = bi,
                  d_neg = d_neg,
                  sigma2e_neg = sigma2e_neg,
                  priors = priors);
    cc_sample(iter,0) = cc;
    //Rcout << "cc" << "\n";
    /////////////////////
    // update mu_a
    ////////////////////
    mu_a = updateMuA(dat = dat,
                     ai = ai,
                     sigma2a = sigma2a,
                     priors = priors);
    mu_a_sample(iter,0) = mu_a;
    //Rcout << "mua" << "\n";

    /////////////////////
    // update mu_b
    ////////////////////
    mu_b = updateMuB(dat = dat,
                     bi = bi,
                     sigma2b = sigma2b,
                     priors = priors);
    mu_b_sample(iter,0) = mu_b;
    //Rcout << "mub" << "\n";
    /////////////////////
    // update sigma2a
    ////////////////////
    sigma2a = updateSigma2a(dat = dat,
                            ai = ai,
                            mu_a = mu_a,
                            priors = priors);
    sigma2a_sample(iter, 0) = sigma2a;
    //Rcout << "s2a" << "\n";
    /////////////////////
    // update sigma2b
    ////////////////////
    sigma2b = updateSigma2b(dat = dat,
                            bi = bi,
                            mu_b = mu_b,
                            priors = priors);
    sigma2b_sample(iter, 0) = sigma2b;
    //Rcout << "s2b" << "\n";
    /////////////////////
    // update sigma2d_neg
    ////////////////////
    sigma2d_neg = updateSigma2Dneg(dat = dat,
                                   d_neg = d_neg,
                                   priors = priors);
    sigma2d_neg_sample(iter, 0) = sigma2d_neg;
    //Rcout << "s2dn" << "\n";
    /////////////////////
    // update sigma2d
    ////////////////////
    sigma2d = updateSigma2D(dat = dat,
                            d_pos = d_pos,
                            d_hk = d_hk,
                            d_reg = d_reg,
                            priors = priors);
    sigma2d_sample(iter,0) = sigma2d;
    //Rcout << "s2d" << "\n";
    /////////////////////
    // update sigma2e_neg
    ////////////////////
    sigma2e_neg = updateSigma2Eneg(dat = dat,
                                   ai = ai,
                                   bi = bi,
                                   cc = cc,
                                   d_neg = d_neg,
                                   priors = priors);
    sigma2e_neg_sample(iter, 0) = sigma2e_neg;
    //Rcout << "s2eneg" << "\n";
    /////////////////////
    // update sigma2e
    ////////////////////
    sigma2e = updateSigma2E(dat = dat,
                            ai = ai,
                            bi = bi,
                            pos_conc = pos_conc,
                            phi = phi,
                            kappa_hk = kappa_hk,
                            kappa_reg = kappa_reg,
                            d_pos = d_pos,
                            d_hk = d_hk,
                            d_reg = d_reg,
                            priors = priors);
    sigma2e_sample(iter, 0) = sigma2e;
    //Rcout << "s2e" << "\n";
    /////////////////////
    // update sigma2kappa_hk
    ////////////////////
    sigma2kappa_hk = updateSigma2KappaHk(dat = dat,
                                         kappa_hk = kappa_hk,
                                         lambda_hk = lambda_hk,
                                         priors = priors);
    sigma2kappa_hk_sample(iter,0) = sigma2kappa_hk;
    //Rcout << "s2khk" << "\n";

    /////////////////////
    // update sigma2kappa_reg
    ////////////////////
    sigma2kappa_reg = updateSigma2KappaReg(dat = dat,
                                         kappa_reg = kappa_reg,
                                         lambda_reg = lambda_reg,
                                         priors = priors);
    sigma2kappa_reg_sample(iter,0) = sigma2kappa_reg;
    //Rcout << "s2kreg" << "\n";
  }
  
  List location_par = List::create(Named("ai") = ai_sample,
                                   Named("bi") = bi_sample,
                                   Named("phi") = phi_sample,
                                   Named("kappa_hk") = kappa_hk_sample,
                                   Named("kappa_reg") = kappa_reg_sample,
                                   Named("lambda_hk") = lambda_hk_sample,
                                   Named("lambda_reg") = lambda_reg_sample,
                                   Named("d_neg") = d_neg_sample,
                                   Named("d_pos") = d_pos_sample,
                                   Named("d_hk") = d_hk_sample,
                                   Named("d_reg") = d_reg_sample,
                                   Named("cc") = cc_sample,
                                   Named("mu_a") = mu_a_sample,
                                   Named("mu_b") = mu_b_sample);
  
  List variance_par = List::create(Named("sigma2a") = sigma2a_sample,
                                   Named("sigma2b") = sigma2b_sample,
                                   Named("sigma2d_neg") = sigma2d_neg_sample,
                                   Named("sigma2d") = sigma2d_sample,
                                   Named("sigma2e_neg") = sigma2e_neg_sample,
                                   Named("sigma2e") = sigma2e_sample,
                                   Named("sigma2kappa_hk") = sigma2kappa_hk_sample,
                                   Named("sigma2kappa_reg") = sigma2kappa_reg_sample);
  
  List result = List::create(Named("location") = location_par,
                             Named("variance") = variance_par);
  return result;
}
