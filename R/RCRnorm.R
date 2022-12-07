#' @title 
#' RCRnorm
#' @export

RCRnorm <- function(dat, 
                    name_of_negative_control_gene_df = "neg_dat",
                    name_of_positive_control_gene_df = "pos_dat",
                    name_of_housekeeping_gene_df = "hk_dat",
                    name_of_regular_gene_df = "reg_dat",
                    pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125)),
                    prior_list = list(m_ab = 9, mm = 3),
                    random_init = FALSE,
                    fast_method = FALSE,
                    n_iter = 8000,
                    burn_in = 5000,
                    n_chain = 4,
                    seed = 42,
                    transform_fun = log10Plus1,
                    ...)
{
  # We first perform a sanity check on the data
  # then we transform the data 
  dat <- dataPreparation(dat, 
                         name_of_positive_control_gene_df,
                         name_of_negative_control_gene_df,
                         name_of_housekeeping_gene_df,
                         name_of_regular_gene_df,
                         transform_fun,
                         ...)
  if(random_init)
  {
    fast_method <- FALSE
  }
  set.seed(seed)
  posterior_sample_list = vector(mode = "list", length = n_chain)
  for(chain in 1 : n_chain)
  {
    # Use ad hoc method to get the model estimates
    initial_values <- getAdHocPointEstimates(dat, pos_conc, random_init)
    # If use fast method, then this is it...
    if(fast_method)
    {
      return(initial_values)
    }
    # Otherwise, we starting setting priors 
    priors <- setPriors(prior_list, initial_values)
    
    # Now we start Gibbs sampling 
    cat("Chain", chain, "\n")
    posterior_sample_list[[chain]] <- sampler(lapply(dat, as.matrix), 
                                              priors,
                                              initial_values,
                                              pos_conc,
                                              n_iter)
  }
  
  return(posterior_sample_list)
}