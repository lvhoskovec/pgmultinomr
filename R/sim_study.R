#' Run a simulation study using pgmultinom
#'
#' @param simnum simulation number 
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in iterations
#' @param n sample size
#' @param miss_prob proportion of cases with missing outcome data 
#' @param allmiss logical; if TRUE then all outcomes are missing for any case with missing outcome data, default is FALSE 
#' @param null_scenario logical; if TRUE then exposures and covariates have no an effect on outcome, default is FALSE
#' @param equal_probs logical; if TRUE the simulation setting is equal probabilities, if FALSE the simulation setting is data probabilities
#'
#' @importFrom mvnfast rmvn
#'
#' @return a list with elements 
#' @export
#'
#'

sim_study = function(simnum = 1, niter=1000, nburn = 500, n=1000, miss_prob = 0,
                     allmiss = FALSE, null_scenario = FALSE, equal_probs = FALSE){
  
  df_all = NULL
  df_cc = NULL
  print(paste("simnum =", simnum))
  set.seed(21 + 129*simnum)

  ############################################
  ### simulate exposure and covariate data ###
  ############################################

  # exposures
  p = 3
  covX = matrix(c(1, -0.13, 0.02, -0.13, 1, -0.86, 0.02, -0.86, 1), 3, 3)
  x = rmvn(n, mu = rep(0,p), sigma = covX) 

  # covariates 
  q = 5
  w_cov = rmvn(n, mu = rep(0,q), sigma = diag(q))
  w = cbind(1, w_cov)
  q = ncol(w)
  
  #############################
  ### simulate outcome data ###
  #############################

  K = 6
  
  ## stick-breaking simulation ##
  if(!null_scenario){
    beta_true = matrix(rnorm(p*(K-1),0,1), ncol = K-1, nrow = p); beta_true # p by K-1
    gamma_true = matrix(rnorm(q*(K-1),0,1), ncol = K-1, nrow = q); gamma_true # q by K-1
    if(equal_probs) gamma_true[1,] = c(-2.8,-2.5,-2.2,-1.3,-0.5) else gamma_true[1,] = c(1.8,0.5,0,0,0) 
  }else{
    beta_true = matrix(0, ncol = K-1, nrow = p); beta_true # p by K-1
    gamma_true = matrix(0, ncol = K-1, nrow = q); gamma_true # q by K-1
    if(equal_probs) gamma_true[1,] = c(-1.8,-1.5,-1,-.5,-.3) else gamma_true[1,] = c(1.2,0.8,0.6,0.2,0.1)
  }

  xbetaK_true = crossprod(t(x),beta_true) # n by K-1
  wgammaK_true = crossprod(t(w),gamma_true) # n by K-1
  psi_true = xbetaK_true + wgammaK_true

  # make pi
  pitildek = exp(psi_true)/(1+exp(psi_true));colMeans(pitildek)
  piik = matrix(0, n, K)
  piik[,1] = pitildek[,1]
  for(k in 2:(K-1)){
    piik[,k] = pitildek[,k]*(1-rowSums(matrix(piik[,(1:(k-1))],nrow=n,ncol = k-1)))
  }
  piik[,K] = 1-rowSums(piik)

  # simulate multinomial data
  # make sure each category has at least one known observation 
  repeat{
    y=t(sapply(1:n, FUN = function(i){
      return(rmultinom(n=1,size=1,prob=piik[i,]))
    }))
    if( !(0 %in% colMeans(y, na.rm = TRUE)) ) break
  }
  
  # save complete data 
  ycomplete = y

  ###################
  ### holdout set ###
  ###################

  if(miss_prob > 0 & !allmiss){
    nmiss = floor(miss_prob*n); nmiss # number of observations missing something
    whichmiss = sort(sample(1:n, nmiss, replace = FALSE)) # which observations are uncertain
    for(i in whichmiss){
      # at least 2 missing for uncertain outcomes 
      which1 = which(ycomplete[i,]==1)
      num_mis = sample(1:(K-1),1)
      extra_mis = sample(which(ycomplete[i,]==0), size = num_mis, replace = FALSE)
      catmiss = sort(c(which1, extra_mis));catmiss
      y[i,catmiss] = NA
    }
  }else if(miss_prob > 0 & allmiss){
    nmiss = floor(miss_prob*n); nmiss # number of observations missing something
    whichmiss = sort(sample(1:n, nmiss, replace = FALSE)) # which observations are uncertain
    y[whichmiss, ] = rep(NA, K)
  }
  
  #################
  ### Fit Model ###
  #################
  
  fit = pgmultinom(niter = niter, priors = NULL, y = y, x = x, w = w,
                   intercept = TRUE)
  
  fit_summary = get_sim_results(fit = fit, ycomplete = ycomplete, beta_true = beta_true, gamma_true = gamma_true)
  
  
  if(miss_prob > 0){
    # complete case fit 
    y_cc = y[-whichmiss,]
    x_cc = x[-whichmiss,]
    w_cc = w[-whichmiss,]

    fit_cc = pgmultinom(niter = niter, priors = NULL, y = y_cc, x = x_cc, w = w_cc,
                        intercept = TRUE)
    
    fit_cc_summary = get_sim_results(fit = fit_cc, beta_true = beta_true, gamma_true = gamma_true)
    
    
  }


  ###################
  ### Imputations ###
  ###################
  
  if(miss_prob > 0){
    precision_mean = apply(fit_summary$precision[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
    recall_mean = apply(fit_summary$recall[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
  }else{
    precision_mean = rep(NA, K)
    recall_mean = rep(NA, K)
  }
  
  ########################################
  ### Regression coefficient estimates ###
  ########################################
  
  # RMSE
  rmse_beta = mean(fit_summary$rmse_beta[(nburn+1):niter])
  rmse_gamma = mean(fit_summary$rmse_gamma[(nburn+1):niter])
  
  # bias
  bias_beta = mean(fit_summary$bias_beta[(nburn+1):niter])
  bias_gamma = mean(fit_summary$bias_gamma[(nburn+1):niter]) 
  
  # coverage 
  beta_quants = t(apply(fit$beta.vec[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), quantile(b, 0.975)))
  }))
  wid_beta = mean(beta_quants[,2]-beta_quants[,1]) # width of interval 
  cov_beta = mean(as.vector(beta_true) > beta_quants[,1] & as.vector(beta_true) < beta_quants[,2])
  
  gamma_quants = t(apply(fit$gamma.vec[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), quantile(b, 0.975)))
  }))
  wid_gamma = mean(gamma_quants[,2]-gamma_quants[,1])
  cov_gamma = mean(as.vector(gamma_true) > gamma_quants[,1] & as.vector(gamma_true) < gamma_quants[,2])

  #############################
  ### save table of results ###
  #############################
  
  df1 = data.frame(rbind(c(rmse_beta, bias_beta, wid_beta, cov_beta, 
                          rmse_gamma, bias_gamma, wid_gamma, cov_gamma,
                          precision_mean, recall_mean)))

  colnames(df1) = c("RMSE_beta", "bias_beta", "ci_width_beta", "cov_beta", 
             "RMSE_gamma", "bias_gamma", "ci_width_gamma", "cov_gamma", 
             "precision1", "precision2", "precision3", "precision4", "precision5", "precision6", 
             "recall1","recall2","recall3","recall4","recall5","recall6")
  
  df_all = rbind(df_all, df1)
  
  
  ######################################################
  ### Complete Case regression coefficient estimates ###
  ######################################################

  if(miss_prob > 0){
    
    # RMSE
    cc_rmse_beta = mean(fit_cc_summary$rmse_beta[(nburn+1):niter])
    cc_rmse_gamma = mean(fit_cc_summary$rmse_gamma[(nburn+1):niter])
    
    # bias
    cc_bias_beta = mean(fit_cc_summary$bias_beta[(nburn+1):niter])
    cc_bias_gamma = mean(fit_cc_summary$bias_gamma[(nburn+1):niter]) 
    
    # coverage
    cc_beta_quants = t(apply(fit_cc$beta.vec[(nburn+1):niter,], 2, FUN = function(b){
      return(c(quantile(b, 0.025), quantile(b, 0.975)))
    }))
    cc_wid_beta = mean(cc_beta_quants[,2]-cc_beta_quants[,1]) # width of interval 
    cc_cov_beta = mean(as.vector(beta_true) > cc_beta_quants[,1] & as.vector(beta_true) < cc_beta_quants[,2])
    
    cc_gamma_quants = t(apply(fit_cc$gamma.vec[(nburn+1):niter,], 2, FUN = function(b){
      return(c(quantile(b, 0.025), quantile(b, 0.975)))
    }))
    cc_wid_gamma = mean(cc_gamma_quants[,2]-cc_gamma_quants[,1])
    cc_cov_gamma = mean(as.vector(gamma_true) > cc_gamma_quants[,1] & as.vector(gamma_true) < cc_gamma_quants[,2])
  }
  
  ###############################################
  ### save table of results for complete case ###
  ###############################################
  
  if(miss_prob > 0){
    df2 = data.frame(rbind(c(cc_rmse_beta, cc_bias_beta, cc_wid_beta, cc_cov_beta, 
                             cc_rmse_gamma, cc_bias_gamma, cc_wid_gamma, cc_cov_gamma)))
    
    colnames(df2) = c("RMSE_beta", "bias_beta", "ci_width_beta", "cov_beta", 
                      "RMSE_gamma", "bias_gamma", "ci_width_gamma", "cov_gamma")
    
    df_cc = rbind(df_cc, df2)
  }
  
  l1 = list(results = df_all, complete_case_results = df_cc)
  return(l1)

}








