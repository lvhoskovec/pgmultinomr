#' Get parameter estimation and imputation results from a simulation 
#'
#' @param fit object of type "pgmultinom" 
#' @param ycomplete (optional) complete outcome data 
#' @param beta_true list of true exposure regression coefficients for each category k = 1,...,K-1
#' @param gamma_true (optional) list of true covariate regression coefficients for each category k = 1,...,K-1
#'
#' @return a list with components
#' \itemize{
#'    \item rmse_beta: RMSE for beta 
#'    \item bias_beta: bias for beta 
#'    \item rmse_gamma: RMSE for gamma if gamma_true is given
#'    \item bias_gamma: bias for gamma if gamma_true is given
#'    \item precision: precision for imputations 
#'    \item recall: recall for imputations 
#' }
#'

get_sim_results = function(fit, ycomplete = NULL, beta_true, gamma_true = NULL){
  
  # number of iterations 
  niter = length(fit$y.save)
  
  # fixed parameters
  n = nrow(fit$y_original)
  K = ncol(fit$y_original)
  ymiss = matrix(as.numeric(is.na(fit$y_original)), nrow = n) 
  nmiss = length(which(rowSums(ymiss)>0))
  
  # estimate parameters
  rmse_beta = numeric()
  rmse_gamma = numeric()
  bias_beta = numeric()
  bias_gamma = numeric()
  
  # impute missing data 
  if(!is.null(ycomplete)){
    precision = matrix(NA, niter, K)
    recall = matrix(NA, niter, K)
    uncertain_k = list()
    k_truths = list()
    total_k = list()
    
    for(k in 1:K){
      uncertain_k[[k]] = which(ymiss[,k] == 1) # which ones are missing for each k 
      k_truths[[k]] = ycomplete[,k][uncertain_k[[k]]]
      total_k[[k]] = length(uncertain_k[[k]])
    }
    
  }else{
    precision = NULL
    recall = NULL
  }
  
  
  # get estimates from fit 
  for(s in 1:niter){
    
    # evaluate classification #
    if(!is.null(ycomplete) & nmiss > 0){
      for(k in 1:K){
        k_preds = fit$y.save[[s]][,k][uncertain_k[[k]]] # predicted values for this iteration and class
        # k_preds = y[,k]
        # true positives 
        tp = length(which(k_preds == 1 & k_truths[[k]] == 1)); tp
        # true negatives 
        tn = length(which(k_preds == 0 & k_truths[[k]] == 0)); tn
        # false positives 
        fp = length(which(k_preds == 1 & k_truths[[k]] == 0)); fp
        # false negatives
        fn = length(which(k_preds == 0 & k_truths[[k]] == 1)); fn
        # precision
        if(tp+fp == 0){
          precision[s,k] = NA
        }else{
          precision[s,k] = tp/(tp + fp)
        }
        # recall 
        if(tp + fn == 0){
          recall[s,k] = NA
        }else{
          recall[s,k] = tp/(tp + fn)
        }
      }
    }
    
  
    # parameter estimation 
    if(!is.null(beta_true)){
      
      betak = fit$beta[[s]]
      
      rmse_beta[s] = sqrt(mean((betak-beta_true)^2))
      bias_beta[s] = mean(betak-beta_true)
    }
    
    if(!is.null(gamma_true)){
      
      gammak = fit$gamma[[s]]
      
      rmse_gamma[s] = sqrt(mean((gammak-gamma_true)^2))
      bias_gamma[s] = mean(gammak-gamma_true)
    }
    
  }
  
  list.save = list(rmse_beta = rmse_beta, bias_beta = bias_beta, 
               rmse_gamma = rmse_gamma, bias_gamma = bias_gamma, 
               precision = precision, recall = recall)
  
  return(list.save)
  
}