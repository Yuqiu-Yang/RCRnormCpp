#' @title 
#' Data preparation
#' @description 
#' This function checks if the input is of the correct format
#' and transforms the data using the user supplied function
#' @param dat A list containing data for the 4 probe types: positive control, negative control, housekeeping gene and regular gene.
#'The names for the 4 elements in the list by default are: pos_dat, neg_dat, hk_dat and reg_dat, respectively. For an example
#'of the input data format, please refer to the FFPE_dat included in the dataset.
#'The data for each probe type should be a dataframe with rows being genes and column being patients. The number of columns (patients)
#'should be the same for data of all four probe types. The rows of positive control data should have the same order as the postive control
#'RNA amount vector supplied to the function.
#' @param name_of_
#' @export

dataPreparation <- function(dat, 
                            name_of_positive_control_gene_df,
                            name_of_negative_control_gene_df,
                            name_of_housekeeping_gene_df,
                            name_of_regular_gene_df,
                            transform_fun = log10Plus1,
                            ...)
{
  if((!is.list(dat)) | 
     (length(dat) != 4))
  {
    stop("The input data should be a list of length 4.")
  }
  n_patient <- unique(sapply(dat, function(x) dim(x)[2]))
  if(length(n_patient) != 1)
  {
    stop("Not identical number of columns (patients).")
  }
  
  gene_types <- c(name_of_positive_control_gene_df,
                  name_of_negative_control_gene_df,
                  name_of_housekeeping_gene_df,
                  name_of_regular_gene_df)
  if((!is.character(gene_types)) | 
     (length(gene_types) != 4))
  {
    stop("We expect a string for the name of the dataframe")
  }
  
  temp <- names(dat)
  temp1 <- setdiff(gene_types, temp)
  if(length(temp1) != 0)
  {
    cat("We expect the names of the 4 dataframes to be", gene_types,
        sep = paste0(c(":", rep(",",3)), c("\n")))
    stop(cat("We can not find a match for",temp1,
             sep = paste0(c(":", rep(",",length(temp1))), c("\n"))))
  }
  
  temp <- c("pos", "neg", "hk", "reg")
  for(i in 1 : 4)
  {
    eval(parse(text = paste(temp[i],"_dat <- transform_fun(dat[[gene_types[i]]],...)",
                            sep = "")))
  }
  return(list(pos_dat = pos_dat,
              neg_dat = neg_dat,
              hk_dat = hk_dat,
              reg_dat = reg_dat))
}
