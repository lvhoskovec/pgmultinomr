#' Fit Bayesian logistic regression model with polya-gamma data augmentation
#'
#' @param niter number of iterations
#' @param priors list of priors
#' @param y matrix of outcome data, with possible missing observations denoted by NA
#' @param ycomplete matrix of complete outcome data, for validation 
#' @param x matrix of exposure data 
#' @param w (optional) matrix of covariate data 
#' @param intercept logical; include an intercept?
#' @param beta_true (optional) list of true betas (exposure coefficients) for simulation study 
#' @param gamma_true (optional) list of true gammas (covariate coefficients) for simulation study 
#' 
#'
#' @return list of parameter estimates
#' @export
#'
pglogistic <- function(niter, priors=NULL, 
                       y, ycomplete=NULL, x, w = NULL, intercept = TRUE,
                       beta_true = NULL, gamma_true = NULL){
  
  n = length(y); n
  q = dim(x)[2]; q
  if(intercept & is.null(w)){
    w = matrix(1, ncol = n)
  }else if(intercept & mean(w[,1]!=1)){
    w = cbind(1, w) # add a column of 1's for the intercept 
  } 
  d = ncol(w); d

  # missing outcome data 
  whichmiss = which(is.na(y))
  nmiss = length(whichmiss)

  ##############
  ### priors ###
  ##############
  
  if(missing(priors)) priors <- NULL
  
  # exposure priors
  if(is.null(priors$betamean)) priors$betamean = rep(0,q)
  if(is.null(priors$betavar)) priors$betavar = diag(q)
  priors$betavarinverse = chol2inv(chol(priors$betavar))
  
  # covariate priors 
  if(!is.null(w)) {
    if(is.null(priors$gammamean)) priors$gammamean = rep(0,d)
    if(is.null(priors$gammavar)) priors$gammavar = diag(d)
    priors$gammavarinverse = chol2inv(chol(priors$gammavar))
  }
  
  
  ##################
  ### imputation ###
  ##################
  
  # impute starting values uniformly
  for(i in 1:n){
    if(is.na(y[i])){
      y[i] = rbinom(1,1,0.5)
    }
  }

  
  ########################
  ### fixed parameters ###
  ########################

  ni = 1
  ytrue = ycomplete[whichmiss]

  #######################
  ### starting values ###
  #######################
  
  beta = rmvn(1, mu = priors$betamean, sigma = priors$betavar)
  gamma = rmvn(1, mu = priors$gammamean, sigma = priors$gammavar)

  xbeta = tcrossprod(x,beta)
  if(!is.null(w)) wgamma = tcrossprod(w, gamma) else wgamma = matrix(0, n, 1)
  psi = xbeta + wgamma
  kappa = y - 1/2
  
  ##############################
  ### space to store results ###
  ##############################
  
  beta.save = matrix(NA, niter, q) # exposure coefficients
  gamma.save = matrix(NA, niter, d) # covariate coefficients 
  y.save = matrix(NA, niter, n) # save imputations 
  rmse_beta = numeric()
  rmse_gamma = numeric()
  bias_beta = numeric()
  bias_gamma = numeric()
  
  if(!is.null(ycomplete)){
    # classification 
    accuracy = numeric(niter)
    precision = numeric(niter)
    recall = numeric(niter)
    specificity = numeric(niter)
    f1score = numeric(niter)
  }else{
    accuracy = NULL
    precision = NULL
    recall = NULL
    specificity = NULL
    f1score = NULL
  }
  
  
  ##################
  ### START MCMC ###
  ##################
  
  for(s in 1:niter){
    
    if(s%%100==0){
      print(s)
    }

    #########################
    ### update parameters ###
    #########################
    
    # update omega
    omega = rcpp_pgdraw(rep(1,n), psi)
    
    # update beta 
    params = update_beta_logistic(omega = omega, kappa = kappa, x = x, wgamma = wgamma,
                                  betaVarInverse = priors$betavarinverse, 
                                  betaMean = priors$betamean)
    beta = rmvn(1, mu = params[[1]], sigma = params[[2]])
    
    # update gamma
    params_gamma = update_gamma_logistic(omega = omega, kappa = kappa, w = w, xbeta = xbeta, 
                                         gammaVarInverse = priors$gammavarinverse, 
                                         gammaMean = priors$gammamean)  
    
    gamma = rmvn(1, mu = params_gamma[[1]], sigma = params_gamma[[2]])
    
    # update xbeta and wgamma
    xbeta = tcrossprod(x,beta)
    if(!is.null(w)) wgamma = tcrossprod(w, gamma) else wgamma = matrix(0, n, 1)
    psi = xbeta + wgamma
    
    ###################################
    ### impute missing outcome data ###
    ###################################
    
    if(nmiss > 0){
      pii = exp(psi)/(1+exp(psi))
      for(i in whichmiss){
        y[i] = rbinom(n = 1, size = 1, prob = pii[i])
      }
    }
    
    kappa = y - 1/2
    
    ###############################
    ### evaluate classification ###
    ###############################
    
    if(!is.null(ycomplete) & nmiss > 0){
      ypreds = y[whichmiss] 

      tp = length(which(ypreds == 1 & ytrue == 1)); tp
      tn = length(which(ypreds == 0 & ytrue == 0)); tn
      fp = length(which(ypreds == 1 & ytrue == 0)); fp
      fn = length(which(ypreds == 0 & ytrue == 1)); fn
      
      accuracy[s] = (tp+tn)/nmiss
      precision[s] = tp/(tp+fp)
      recall[s] = tp/(tp + fn)
      specificity[s] = tn/(tn + fp)
      f1score[s] = 2*(precision[s] * recall[s])/(precision[s] + recall[s])
    }
    
  
    ######################################
    ### RMSE and bias for coefficients ###
    ######################################
    
    if(!is.null(beta_true)){
      rmse_beta[s] = sqrt(mean((beta-beta_true)^2))
      bias_beta[s] = mean(beta-beta_true)
    }
    
    if(!is.null(gamma_true)){
      rmse_gamma[s] = sqrt(mean((gamma-gamma_true)^2))
      bias_gamma[s] = mean(gamma-gamma_true)
    }
    
  
    ######################
    ### save estimates ###
    ######################
    
    y.save[s,] = y
    beta.save[s,] = beta
    gamma.save[s,] = gamma
    

  }
  
  list.save = list(beta = beta.save, gamma = gamma.save, 
                   rmse_beta = rmse_beta, bias_beta = bias_beta, 
                   rmse_gamma = rmse_gamma, bias_gamma = bias_gamma, 
                   y.save = y.save,
                   accuracy = accuracy, precision = precision, recall = recall, 
                   specificity = specificity, f1score = f1score)
  
  class(list.save) = "pglogistic"
  
  return(list.save)
}