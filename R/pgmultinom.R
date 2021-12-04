#' Bayesian categorical regression model with polya-gamma data augmentation
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
#' @importFrom mvnfast rmvn 
#' @importFrom stats rmultinom
#'
#' @return a list with components
#' \itemize{
#'    \item beta: list with components the matrix of estimated exposure regression coefficients at each iteration
#'    \item beta.vec: list with components vectorized estimated exposure regression coefficients at each iteration
#'    \item gamma: list with components the matrix of estimated covariate regression coefficients at each iteration
#'    \item gamma.vec: list with components the vectorized estimated covariate regression coefficients at each iteration
#'    \item rmse_beta: RMSE for beta if beta_true is given
#'    \item bias_beta: bias for beta if beta_true is given
#'    \item rmse_gamma: RMSE for gamma if gamma_true is given
#'    \item bias_gamma: bias for gamma if gamma_true is given
#'    \item y.save: list with components the matrix of outcomes with imputations
#'    \item precision: precision for each category, if ycomplete is given
#'    \item recall: recall for each category, if ycomplete is given 
#' 
#' 
#' }
#' @export
#'
pgmultinom <- function(niter, priors=NULL, 
                       y, ycomplete=NULL, x, w = NULL, intercept = TRUE,
                       beta_true = NULL, gamma_true = NULL){
  
  
  K = dim(y)[2]; K # number of categories
  n = dim(y)[1]; n # sample size
  q = dim(x)[2]; q # number of exposures 
  
  if(intercept & is.null(w)){
    w = matrix(1, ncol = n)
  }else if(intercept & mean(w[,1]!=1)){
    w = cbind(1, w) # add a column of 1's for the intercept 
  } 
  d = dim(w)[2]; d # number of covariates 
  
  # missing outcome data 
  ymiss = matrix(as.numeric(is.na(y)), nrow = n) 
  nmiss = length(which(rowSums(ymiss)>0))
  whichmiss = which(rowSums(ymiss)>0)
  
  ##############
  ### priors ###
  ##############
  
  if(missing(priors)) priors <- NULL
  
  # exposure priors
  if(is.null(priors$betamean)) priors$betamean = lapply(1:(K-1), FUN = function(k) rep(0,q))
  if(is.null(priors$betavar)) priors$betavar = lapply(1:(K-1), FUN = function(k) diag(q))
  priors$betavarinverse = lapply(1:(K-1), FUN = function(k) chol2inv(chol(priors$betavar[[k]])))
  
  # covariate priors 
  if(!is.null(w)) {
    if(is.null(priors$gammamean)) priors$gammamean = lapply(1:(K-1), FUN = function(k) rep(0,d))
    if(is.null(priors$gammavar)) priors$gammavar = lapply(1:(K-1), FUN = function(k) diag(d))
    priors$gammavarinverse = lapply(1:(K-1), FUN = function(k) chol2inv(chol(priors$gammavar[[k]])))
  }


  ##################
  ### imputation ###
  ##################
  
  # impute starting values uniformly
  for(i in 1:n){
    if(sum(ymiss[i,])>0){
      k_select = which(ymiss[i,] == 1) # which categories are missing 
      y[i,k_select] = rmultinom(1, size = 1, prob = rep(1/length(k_select), length(k_select)))
    }
  }

  ########################
  ### fixed parameters ###
  ########################
  
  nik = matrix(0,n,K-1)
  nik[,1] = 1
  for(k in 2:(K-1)){
    ypartsum = matrix(y[,(1:(k-1))],nrow=n, ncol=k-1)
    nik[,k] = 1 - rowSums(ypartsum)
  }

  whichik = lapply(1:n, FUN = function(i){
    which(nik[i,]==1)
  })

  # make the pieces to update beta 
  the_xs = list() # submatrix of X containing the rows i such that nik = 1
  the_ws = list() # ditto for w 
  the_ys = list()
  omega = list()

  for(k in 1:(K-1)){
    sub = which(nik[,k]==1)
    the_xs[[k]] = matrix(x[sub,], ncol = q)
    if(!is.null(w)) the_ws[[k]] = matrix(w[sub,], ncol = d)
    the_ys[[k]] = y[sub,k]-(1/2)
  }
  
  #######################
  ### starting values ###
  #######################
  
  betak = matrix(rnorm(q*(K-1)), nrow = q, ncol = K-1)
  if(!is.null(w)) gammak = matrix(rnorm(d*(K-1)), nrow = d, ncol = K-1) 

  ##############################
  ### space to store results ###
  ##############################
  
  beta.save = list()# exposure coefficients
  gamma.save = list() # covariate coefficients 
  y.save = list() # save imputations 
  rmse_beta = numeric()
  rmse_gamma = numeric()
  bias_beta = numeric()
  bias_gamma = numeric()
  beta.vec.save = matrix(NA, niter, q*(K-1))
  gamma.vec.save = matrix(NA, niter, d*(K-1))

  if(!is.null(ycomplete)){
    # classification 
    accuracy = matrix(NA, niter, K)
    precision = matrix(NA, niter, K)
    recall = matrix(NA, niter, K)
    specificity = matrix(NA, niter, K)
    f1score = matrix(NA, niter, K)
    uncertain_k = list()
    k_truths = list()
    total_k = list()
    
    for(k in 1:K){
      uncertain_k[[k]] = which(ymiss[,k] == 1) # which ones are missing for each k 
      k_truths[[k]] = ycomplete[,k][uncertain_k[[k]]]
      total_k[[k]] = length(uncertain_k[[k]])
    }
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

    #################################################
    ### parameters that can change when we impute ###
    #################################################
    
    nik = matrix(0,n,K-1)
    nik[,1] = 1
    for(k in 2:(K-1)){
      ypartsum = matrix(y[,(1:(k-1))],nrow=n)
      nik[,k] = 1 - rowSums(ypartsum)
    }
    
    whichik = lapply(1:n, FUN = function(i){
      which(nik[i,]==1)
    })
    
    # make the pieces to update beta 
    the_xs = list() # submatrix of X containing the rows i such that nik = 1
    the_ws = list()
    the_ys = list()
    omega = list()
    # will change when we impute 
    for(k in 1:(K-1)){
      sub = which(nik[,k]==1)
      the_xs[[k]] = matrix(x[sub,], ncol = q)
      if(!is.null(w)) the_ws[[k]] = matrix(w[sub,], ncol = d)
      the_ys[[k]] = y[sub,k]-(1/2)
    }
    
    #########################################
    ### calculate x times beta for each k ###
    #########################################
    
    xbetak = x %*% betak
    if(!is.null(w)) zgammak = w %*% gammak else zgammak = matrix(0, n, K-1)
    
    ######################################
    ### update omega, the PG variables ###
    ######################################
    
    omega = update_omega(n=n, K=K-1, whichik = whichik, 
                                xbetak=xbetak, zgammak = zgammak)
    
    the_omegas = lapply(1:(K-1), FUN= function(k) {
      return(omega[which(omega[,k]!=0),k])
    })
    
    
    ####################
    ### update betak ###
    ####################
    
    # parallelizing the sampling of beta_k 
    params = lapply(1:(K-1), FUN = function(k){  
      if(!is.null(w)){
        wgam_for_beta = the_ws[[k]] %*% gammak[,k]
      }else {
        wgam_for_beta = matrix(0, nrow = length(the_ys[[k]]))
      }
      
      update_beta(k = k-1, ok = the_omegas[[k]], yk = the_ys[[k]], xk = the_xs[[k]], 
                  wgam = wgam_for_beta,
                  betaVarInverse = priors$betavarinverse, betaMean = priors$betamean)
      })
      
    for(k in 1:(K-1)){
      betak[,k] = rmvn(1, mu = params[[k]][[1]], sigma = params[[k]][[2]])
    }


    #####################
    ### update gammaK ###
    #####################
    
    if(!is.null(w)){
      
      params_cov = lapply(1:(K-1), FUN = function(k){
        
        xbet_for_gamma = the_xs[[k]]%*%betak[,k]
        
        update_gamma(k = k-1, ok = the_omegas[[k]], yk = the_ys[[k]], wk = the_ws[[k]], 
                     xbet = xbet_for_gamma,
                     gammaVarInverse = priors$gammavarinverse, gammaMean = priors$gammamean)
      })
      
      for(k in 1:(K-1)){
        gammak[,k] = rmvn(1, mu = params_cov[[k]][[1]], sigma = params_cov[[k]][[2]])
      }
    }

    ###################################
    ### impute missing outcome data ###
    ###################################
    
    if(nmiss > 0){
      psi = xbetak + zgammak
      pitildek = exp(psi)/(1+exp(psi))
      piik = matrix(0, n, K)
      piik[,1] = pitildek[,1]
      for(k in 2:(K-1)){
        piik[,k] = pitildek[,k]*(1-rowSums(matrix(piik[,(1:(k-1))],nrow=n,ncol = k-1)))
      }
      piik[,K] = 1-rowSums(piik) # multinomial probabilities, Kth category
      
      for(i in whichmiss){
        k_select = which(ymiss[i,] == 1) # which categories are missing
        y[i,k_select] = rmultinom(1, size = 1, prob = piik[i,][k_select])
      }
    }

    ###############################
    ### evaluate classification ###
    ###############################

    if(!is.null(ycomplete) & nmiss > 0){
      for(k in 1:K){
        k_preds = y[,k][uncertain_k[[k]]] # predicted values for this iteration and class
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

    ######################################
    ### RMSE and bias for coefficients ###
    ######################################
    
    if(!is.null(beta_true)){
      rmse_beta[s] = sqrt(mean((betak-beta_true)^2))
      bias_beta[s] = mean(betak-beta_true)
    }

    if(!is.null(gamma_true)){
      rmse_gamma[s] = sqrt(mean((gammak-gamma_true)^2))
      bias_gamma[s] = mean(gammak-gamma_true)
    }
    
    
    ################################
    ### vectorize beta and gamma ###
    ################################
    
    if(!is.null(w)) {
      gamma.save[[s]] = gammak
      gamma.vec.save[s,] = as.vector(gammak)
    }
    beta.save[[s]] = betak
    beta.vec.save[s,] = as.vector(betak)
    
    ##############
    ### save y ###
    ##############
    y.save[[s]] = y
    
  }
  
  list.save = list(beta = beta.save, beta.vec = beta.vec.save, 
                   gamma = gamma.save, gamma.vec = gamma.vec.save,
                   rmse_beta = rmse_beta, bias_beta = bias_beta, 
                   rmse_gamma = rmse_gamma, bias_gamma = bias_gamma, 
                   y.save = y.save,
                   precision = precision, recall = recall)
  
  class(list.save) = "pgmultinom"
  
  return(list.save)
}