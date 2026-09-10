# k-fold cross validation function for SSN models
# 
# ssnfitmodel: a fitted SSN model to use for k-fold cross validation. [REQUIRED
#              INPUT FOR FUNCTION]
# 
# kfoldn:      the number of folds/data chunks to use for splitting the data.
#              5- and 10-fold cross validation are typical, but values can vary.
#              Values can only go to 99 due to 2-digit padding in kfold names. 
#              [REQUIRED INPUT FOR FUNCTION]
# 
# obsID:       an ID field used to track observation sites within an SSN object.
#              Using this column allows the function to match observations 
#              withheld for kfold predictions to their original observed values
#              for comparison later. Default is "pid" since all SSN object have
#              this field, but alternative point ID tracking fields can be 
#              specified when desired. [OPTIONAL INPUT FOR FUNCTION]
# 
# clustern:    The number of processors to use in parallel processing of SSN 
#              model fits and predictions for each withheld k-fold (default = 2).
#              If unsure how many processors to choose, set clustern = NA and
#              the function will detect and select a maximum number of cores to 
#              use via: `clustern <- (parallel::detectCores() - 2)`
#              This maximum number of coures may not be ideal if kfoldn is 
#              smaller than the value selected by 'detectCores() - 2'. The
#              function will still run, but there were be idle parallel cores
#              spun up that will be unused. This could slow down the computer.
#              [OPTIONAL INPUT FOR FUNCTION]
# 
# requires:    dplyr, parallel, pbapply, purrr, rsample, sf, SSN2, stats, 
#              stringr, tibble, tidyr


SSNkfoldcv <- 
  function(ssnfitmodel, kfoldn = 5, clustern = 2, obsID = "pid"){
    
    # grab ssn object from fitted SSN model object
    ssnobj <- ssnfitmodel$ssn.object
    
    # create a list of kfold names that are padded with zeros
    kfold_names <- paste0("k", stringr::str_pad(1:kfoldn, width=2, 
                                                side="left", pad="0") )
    
    # response variable (as character string) to set to NA in k-fold cv chunks
    ssn_responsevar <- all.vars(ssnfitmodel$call$formula)[1]
    
    # extract the observation site data frame
    ssn_obs_df00 <- SSN2::ssn_get_data(ssnobj, "obs")
    
    # split the data frame into k-folds as a named list: k01, k02, k03, ...kn
    # with equal numbers of rows or nearly so if not divisible by specified k
    set.seed(12); kfolds_ls <- 
      split(x = sample(ssn_obs_df00[[obsID]], nrow(ssn_obs_df00), replace = FALSE), 
            f = as.factor(kfold_names) )
    
    # Create a data frame key to assigned a k-fold group to each observation
    kfolds_df <- 
      tidyr::unnest(tibble::enframe(kfolds_ls, name = "kfold", value = obsID), 
                    cols = dplyr::all_of(obsID)) |> data.frame()
    
    # Join k-fold groups to ssn observation data frame by the obsID column
    ssn_obs_df01 <- dplyr::left_join(ssn_obs_df00, kfolds_df, by = obsID)
    
    # Function to set response variable to NA for the assigned k-fold (k) group
    # k: is the name of the k-fold provided for tracking (e.g., k01, k02, 
    # k03, ..., kn)
    kfold_na_set_fnx <- function(k) {
      dfkn <- ssn_obs_df01
      dfkn[dfkn$kfold == k, ssn_responsevar] <- NA
      return(dfkn)
    }
    
    # Apply k-fold NA assignment function above and create list of data frames 
    # that hold the model fitting data frames for each k-fold subset
    ssn_obs_dfkn_ls <- lapply(kfold_names, kfold_na_set_fnx)
    names(ssn_obs_dfkn_ls) <- kfold_names  # rename list elements
    
    
    # Function to fit each k-fold withheld data set and make predictions for 
    # observations set to NA. Output from this function is a data frame with 
    # the predictions for each withheld k-fold data set to NA and its ObsID. 
    kfold_ssn_fit_fnx <- function(list_of_dfs) {
      # put the k-fold data frame into the ssn object for updated split/k-fold
      ssn_obj_k <-
        SSN2::ssn_put_data(data = list_of_dfs, x = ssnobj, name = "obs")
      
      # refit the model using the updated ssn object with NAs added for k-folds
      ssnfit <- stats::update(ssnfitmodel, ssn.object = ssn_obj_k)
      
      # make predictions for the missing NA values in the observed data set
      # then select just the k-fold ID, obsID, and fitted/predicted value columns
      ssnpred <- SSN2::augment(ssnfit, newdata = ".missing") |>
        dplyr::select(kfold, dplyr::contains(obsID), .fitted) |>
        sf::st_drop_geometry()
      
      return(ssnpred)
    }
    
    
    # establish parallel cluster to use parallel processing for k-fold 
    # fitting and predictions
    if(is.na(clustern)) {clustern <- (parallel::detectCores() - 2)}
    clust <- parallel::makeCluster(clustern)
    parallel::clusterExport(clust, varlist = c("obsID"), envir = environment())
    parallel::clusterEvalQ(clust, library("SSN2"))
    
    # Apply k-fold fit/prediction function above to each folded data set and 
    # reduce the list of data frames into a single data frame. Then, join to the
    # original data to match observed data with their predictions from the
    # k-fold predictions.
    kfold_preds <- pbapply::pblapply(ssn_obs_dfkn_ls, 
                                     kfold_ssn_fit_fnx,
                                     cl = clust) |>
      purrr::reduce(.f = dplyr::bind_rows) |>
      dplyr::left_join(dplyr::select(ssn_obs_df00,
                                     dplyr::contains(c(obsID,ssn_responsevar))) |>
                         sf::st_drop_geometry(), by = obsID)
    
    # cancel the cluster 
    parallel::stopCluster(clust)
    
    
    # calculate loss function statistics for each fold and overall (entire 
    # data set average across all kfolds). This one-row data frame is the final 
    # output of the function.
    lossfnx_df <- kfold_preds |>
      dplyr::bind_rows((dplyr::mutate(kfold_preds, kfold = "Overall")) ) |>
      dplyr::mutate(error_val = .fitted - .data[[ssn_responsevar]]) |>
      dplyr::mutate(errsq = error_val^2) |>
      dplyr::group_by(kfold) |>
      dplyr::summarise(bias = mean(error_val, na.rm = TRUE),
                       cor2 = cor(.fitted, .data[[ssn_responsevar]])^2,
                       MSPE = mean(errsq, na.rm = TRUE) ) |>
      dplyr::mutate(RMSPE = sqrt(MSPE)) 
    
    return(lossfnx_df)
    
  }
