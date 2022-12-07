#' @title 
#' Get ad hoc point estimates
#' @description 
#' This function uses ANOVA to give quick estimate of all the model parameters 
#' @export

getAdHocPointEstimates <- function(dat,
                             pos_conc, random = F)
{
  pos_dat <- dat$pos_dat
  neg_dat <- dat$neg_dat
  hk_dat <- dat$hk_dat
  reg_dat <- dat$reg_dat
  #number of each class of genes.
  n_hk = dim(hk_dat)[1]
  n_reg = dim(reg_dat)[1]
  n_neg = dim(neg_dat)[1]
  n_pos = dim(pos_dat)[1]
  #number of patients or samples
  n_patient = dim(pos_dat)[2]
  # We first randomly initialize 
  ai <- stats::rnorm(n_patient, 2.5, .1)
  bi <- stats::rnorm(n_patient, .9, .1)
  
  mu_a <- stats::runif(1, 1, 4)
  mu_b <- stats::runif(1, 0, 2)
  sigma2a <- stats::runif(1, 0, .01)
  sigma2b <- stats::runif(1, 0, .01)
  
  cc <- stats::runif(1, -6, -1)
  
  lambda_hk <- stats::rnorm(n_hk, 0, 1)
  lambda_reg <- stats::rnorm(n_reg, 0, 1)
  
  phi <- stats::rnorm(n_patient, 0, 2)
  phi[n_patient] <- -sum(phi[1:(n_patient - 1)])
  
  
  kappa_hk <- stats::rnorm(n_hk * n_patient, 0, 1)
  kappa_reg <- stats::rnorm(n_reg * n_patient, 0, 1)
  sigma2kappa_hk <- stats::runif(1, 0, 1)
  sigma2kappa_reg <- stats::runif(1, 0, 1)
  
  d_neg <- stats::rnorm(n_neg, 0, .01)
  d_pos <- stats::rnorm(n_pos, 0, .01)
  d_hk <- stats::rnorm(n_hk, 0, .01)
  d_reg <- stats::rnorm(n_reg, 0, .01)
  
  sigma2d_neg <- stats::runif(1, 0, .1)
  sigma2d <- stats::runif(1, 0, .1)
  
  sigma2e_neg <- stats::runif(1, 0, .1)
  sigma2e <- stats::runif(1, 0, .1)
  
  phi_sd <- 1
  norm_dat <- pos_RNA <- neg_RNA <- hk_RNA <- reg_RNA <- NULL
  
  if(!random)
  {
    #######################################
    # ai  and   bi
    #######################################
    # We will use simple linear regression to estimate 
    # ai and bi for each patient 
    for(i in 1 : n_patient)
    {
      mod <- lm(pos_dat[,i] ~ pos_conc)
      ai[i] <- mod$coefficients[1]
      bi[i] <- mod$coefficients[2]
    }
    
    
    #######################################
    # mu_a, mu_b, sigma2a, and sigma2b
    #######################################
    # mu_a, mu_b, sigma2a, and sigma2b can then 
    # be estimated using mean and var
    mu_a = mean(ai)
    mu_b = mean(bi)
    sigma2a <- var(ai)
    sigma2b <- var(bi)
    
    #######################################
    # X's and c
    #######################################
    # Now that we have the intercepts and slopes
    # we can use this information to estimate 
    # RNA levels for each type of genes 
    #
    # For housekeeping genes and regular genes
    # we subtract the intercepts from the corresponding 
    # obserations and then devide by the slopes 
    hk_RNA <- sweep(sweep(hk_dat, 2, ai, '-'), 2, bi, '/')
    reg_RNA <- sweep(sweep(reg_dat, 2, ai, '-'), 2, bi, '/')
    # For negative control genes, the same procedure 
    # is used. However, since we believe all negative 
    # controls share the same RNA levels, we will 
    # average the results
    cc <- mean(unlist(sweep(sweep(neg_dat, 2, ai, '-'), 2, bi, '/')))
    # We construct a matrix for future use 
    neg_RNA <- matrix(rep(cc, n_neg * n_patient), ncol = n_patient)
    # Since we know the RNA for positive controls 
    # we simply use that information to construct the matrix 
    pos_RNA <- matrix(rep(pos_conc, n_patient), ncol = n_patient)
    
    #######################################
    # lambda_hk   and    lambda_reg
    #######################################
    # We can then utilize the estimated RNA levels 
    # for the housekeeping genes and regular genes 
    # to estimate the gene-specific effects lambda's 
    lambda_hk <- apply(hk_RNA, 1, mean)
    lambda_reg <- apply(reg_RNA, 1, mean)
    
    #######################################
    # phi
    #######################################
    # estimate patient effect phi by two way ANOVA 
    # with patient's regular gene expression level.
    gene <- factor(rep(1:n_reg, n_patient))
    patient <- factor(rep(1:n_patient, each = n_reg))
    mod <- lm(unlist(reg_RNA) ~ patient + gene, 
              contrasts = list(patient = 'contr.sum', gene = 'contr.sum'))
    norm_dat <- matrix(mod$residuals, nrow = n_reg)
    
    phi <- numeric(n_patient)
    phi[1:(n_patient - 1)] <- summary(mod)$coefficients[2:n_patient, 1]
    phi[n_patient] <- -sum(phi)
    
    phi_sd <- summary(mod)$coefficients[2, 2]
    
    #######################################
    # kappa_hk, kappah_reg, sigma2kappa_hk, and sigma2kappa_reg
    #######################################
    # Now we can combine the estimated RNA levels and 
    # the patient-specific effects phi to estimate 
    # kappa's for both housekeeping genes and regular genes 
    estimate_kappa <- sweep(rbind(hk_RNA, reg_RNA), 2, phi, '-')
    estimate_kappa_var <- sweep(sweep(rbind(hk_RNA, reg_RNA), 2, phi, '-'), 1, c(lambda_hk, lambda_reg), '-')
    
    kappa_hk <- as.vector(unlist(estimate_kappa[1:n_hk,]))
    sigma2kappa_hk <- stats::var(as.vector(unlist(estimate_kappa_var[1:n_hk,])))     #0.02497646
    kappa_reg <- as.vector(unlist(estimate_kappa[(1+n_hk):(n_hk+n_reg),]))
    sigma2kappa_reg <- stats::var(as.vector(unlist(estimate_kappa_var[(1+n_hk):(n_hk+n_reg),])) )   #0.1220459
    
    #######################################
    # d's, sigma2d_neg, and sigma2d
    #######################################
    # We can get the residuals from each type of genes
    # using the estimated RNA levels combined with 
    # the intercepts and slopes. This information 
    # can be used to estimate lane effect d's
    d_neg <- neg_dat - sweep(sweep(neg_RNA, 2, bi, '*'), 2, ai, '+')
    d_neg <- apply(d_neg, 1, mean)
    sigma2d_neg <- var(d_neg)
    
    d_pos <- pos_dat - sweep(sweep(pos_RNA, 2, bi, '*'), 2, ai, '+')
    d_pos <- apply(d_pos, 1, mean)
    sigma2d <- var(d_pos)
    
    d_hk <- hk_dat - sweep(sweep(hk_RNA, 2, bi, '*'), 2, ai, '+')
    d_hk <- apply(d_hk, 1, mean)
    
    d_reg <- reg_dat - sweep(sweep(reg_RNA, 2, bi, '*'), 2, ai, '+')
    d_reg <- apply(d_reg, 1, mean)
    
    #######################################
    # sigma2e_neg   and    sigma2e
    #######################################
    sigma2e_neg <- neg_dat - sweep(sweep(neg_RNA, 2, bi, '*'), 2, ai, '+')
    sigma2e_neg <- sweep(sigma2e_neg, 1, d_neg, '-')
    sigma2e_neg <- var(unlist(sigma2e_neg))
    
    sigma2e <- pos_dat - sweep(sweep(pos_RNA, 2, bi, '*'), 2, ai, '+')
    sigma2e <- sweep(sigma2e, 1, d_pos, '-')
    sigma2e <- var(unlist(sigma2e))
  }
  
  return(list(norm_dat = norm_dat,
              ai = ai,
              bi = bi,
              mu_a = mu_a,
              mu_b = mu_b,
              sigma2a = sigma2a,
              sigma2b = sigma2b,
              pos_RNA = pos_RNA,
              neg_RNA = neg_RNA,
              cc = cc,
              hk_RNA = hk_RNA,
              reg_RNA = reg_RNA,
              phi = phi,
              phi_sd = phi_sd,
              kappa_hk = kappa_hk,
              kappa_reg = kappa_reg,
              lambda_hk = lambda_hk,
              lambda_reg = lambda_reg,
              sigma2kappa_hk = sigma2kappa_hk,
              sigma2kappa_reg = sigma2kappa_reg,
              d_pos = d_pos,
              d_neg = d_neg,
              d_hk = d_hk,
              d_reg = d_reg,
              sigma2d = sigma2d,
              sigma2d_neg = sigma2d_neg,
              sigma2e = sigma2e,
              sigma2e_neg = sigma2e_neg))
}
