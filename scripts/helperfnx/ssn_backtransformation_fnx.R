# Back transformation of covariates from an SSN model fit that used standardize 
# independent variables during the model fitting process. Back transformed 
# variables are then represented via their original units.
# Final estimate output table provides the standardized estimates, raw 
# estimates, standard errors for both estimates, t-statistic, and p-values.
# Basically, we add the back transformed estimates and their SE values to the
# output of 'tidy()' from SSN2.
#
# continuous_vars_df: data frame from observation data that holds the continuous
#                     variables that were original standardized via the 'stand()'
#                     function written for this type of SSN analysis. See:
#                     stand <- function(x) {(x-mean(x, na.rm = T))/(2*sd(x, na.rm = T))}
#                     This function takes that same data frame and repeats the 
#                     standardization process to get a back transformation factor
#                     to recalculate the standardized estimate coefficients 
#                     back to a raw-units estimate coefficient for each variable.
#                     
# tidy_out_tibble: This is the output from running SSN2::tidy() on a fitted
#                     SSN model object. 
#                     
# roundval: the number of digits to round values to across the estimate table.
#                     If zeros are encountered in the table, increasing the
#                     roundval number of digits may provide more detailed values
#                     
# message: logical value to suppress the warning about NA values for the 
#                     intercept when back transforming standardized coefficients
#                     to raw unit coefficients.
#                     
#                     
# NOTE: it is normal and expected for the intercept to have 'NAs" for the raw
#       estimate and raw SE values since we don't have back transformation
#       factors to make that back transformation for it.



std_to_raw_estimate_table_fnx <- 
  function(continuous_vars_df, tidy_out_tibble, roundval = 5, message = TRUE) {
    
    # saved continuous variable data frame with raw values/original units for 
    # each variable from ssn object set up
    continuous <- continuous_vars_df  
    
    # grab just the continuous variables and back transform them to raw values
    # since we fit the model using the standardized values
    # generate the standardization factors to divide the estimates by
    backtrans_2sd <- (2*sapply(continuous,sd)) %>% data.frame() %>% 
      rename(backtransfac = ".") %>%
      rownames_to_column(var = "covariates") %>%
      mutate(covariates = paste0("s",covariates))
    
    # grab the covariates from the continuous variables data frame
    # join these to the back transform factor in same order as in estimate table
    backtrans_est_tab <- tidy_out_tibble %>%
      left_join(backtrans_2sd, by = c("term" = "covariates") ) %>%
      rename(std_est = estimate, std_se = std.error, 
             t_stat = statistic, p_val = p.value) %>%
      mutate(raw_est = std_est/backtransfac, raw_se = std_se/backtransfac) %>%
      select(term, starts_with("std_"), starts_with("raw"), t_stat, p_val) %>%
      mutate(across(where(is.numeric), .fns = function(x) {round(x,roundval)} ) )
    # **********************************************************************
    
    if(isTRUE(message)) {
      message(paste0("It is normal and expected for the intercept to have NAs for ",
                     "the raw estimate (raw_est) and raw SE (raw_se) values since ",
                     "we don't have a back transformation factor to make that ",
                     "transformation for the intercept.") 
              )
    }
    
    
    
    return(backtrans_est_tab)
  }
