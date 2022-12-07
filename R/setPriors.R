#' @title 
#' Set priors
#' @export
#' 

setPriors <- function(user_prior_list,
                      initial_values)
{
  m_ab <- user_prior_list[["m_ab"]]
  mm <- user_prior_list[["mm"]]
  
  if(is.null(m_ab) | is.null(mm))
  {
    stop("Need a value for mm and m_ab.")
  }
  
  priors <- list(u = NULL, 
                 v = NULL,
                 mu_a_hat = NULL,
                 mu_b_hat = NULL,
                 sigma2_mu_a = NULL,
                 sigma2_mu_b = NULL,
                 lambda_hk_lb = NULL,
                 lambda_hk_ub = NULL,
                 lambda_reg_lb = NULL,
                 lambda_reg_ub = NULL,
                 phi_lb = NULL,
                 phi_ub = NULL,
                 cc_lb = NULL,
                 cc_ub = NULL)
  
  for(pname in names(priors))
  {
    temp <- user_prior_list[[pname]]
    if(is.null(temp))
    {
      if(pname %in% c("mu_a_hat", "mu_b_hat", 
                      "sigma2_mu_a", "sigma2_mu_b"))
      {
        n_patient <- length(initial_values$ai)
        mu_a = numeric(500)
        mu_b = numeric(500)
        for (i in 1:500)
        {
          mu_a[i] = mean(sample(initial_values$ai, n_patient - 2))
          mu_b[i] = mean(sample(initial_values$bi, n_patient - 2))
        }
        priors[["mu_a_hat"]] = mean(mu_a)
        priors[["mu_b_hat"]] = mean(mu_b)
        priors[["sigma2_mu_a"]] = m_ab * stats::var(mu_a)
        priors[["sigma2_mu_b"]] = m_ab * stats::var(mu_b)
      }else if(pname %in% c("lambda_hk_lb", "lambda_hk_ub",
                      "lambda_reg_lb", "lambda_reg_ub")){
        priors[["lambda_hk_lb"]] = apply(initial_values$hk_RNA, 1, function(x) return(mean(x) - mm * sd(x)))
        priors[["lambda_hk_ub"]] = apply(initial_values$hk_RNA, 1, function(x) return(mean(x) + mm * sd(x)))
        
        priors[["lambda_reg_lb"]] <- apply(initial_values$reg_RNA, 1, function(x) return(mean(x) - mm * sd(x)))
        priors[["lambda_reg_ub"]] <- apply(initial_values$reg_RNA, 1, function(x) return(mean(x) + mm * sd(x)))
      }else if(pname %in% c("phi_lb", "phi_ub")){
        priors[["phi_lb"]] = initial_values$phi - mm * initial_values$phi_sd
        priors[["phi_ub"]] = initial_values$phi + mm * initial_values$phi_sd
      }else if(pname %in% c("cc_lb", "cc_ub")){
        priors[["cc_lb"]] <- -6
        priors[["cc_ub"]] <- -1
      }else if(pname %in% c("u", "v")){
        priors[["u"]] <- 0.01
        priors[["v"]] <- 0.01
      }

    }else{
      priors[[pname]] <- temp
    }
  }

  return(priors)
}