#' Run a simulation study using PG Bayesian multinomial sampler
#'
#' @param simnum simulation number 
#' @param miss_prob proportion of cases with missing outcome data 
#' @param niter number of MCMC iterations
#' @param priors list of priors
#' @param n sample size
#' @param K number of outcome categories 
#' @param covX covariance of exposure data for simulation 
#' @param allmiss logical; if TRUE then all outcomes are missing for any case with missing outcome data, default is FALSE 
#' @param null_scenario logical; if TRUE then exposures and covariates have no an effect on outcome, default is FALSE
#' @param equal_probs logical; if TRUE the simulation setting is equal probabilities, if FALSE the simulation setting is data probabilities
#'
#' @importFrom mvnfast rmvn
#'
#' @return summary of simulation results
#' @export
#'
#'
#'
sim_study = function(simnum = 1, miss_prob, niter=1000, nburn = 500, priors=NULL, n=1000, K=6,
                     covX = diag(3), allmiss = FALSE, null_scenario = FALSE, equal_probs = FALSE){
  
  df_all = NULL
  df_cc = NULL
  # n = 1000 # sample size
  # K = 6 # number of categories
  
  print(paste("simnum =", simnum))
  set.seed(21 + 129*simnum)

  ############################################
  ### simulate exposure and covariate data ###
  ############################################

  # simulate x data with true covariance structure 
  q = dim(covX)[1]
  x = rmvn(n, mu = rep(0,q), sigma = covX) # change the mean to -1 

  # simple covariates
  d = 5
  w_cov = rmvn(n, mu = rep(0,d), sigma = diag(d))
  w = cbind(1, w_cov)
  d = ncol(w)
  
  #############################
  ### simulate outcome data ###
  #############################

  ## stick-breaking simulation ##
  if(!null_scenario){
    beta_true = matrix(rnorm(q*(K-1),0,1), ncol = K-1, nrow = q); beta_true # q by K-1
    gamma_true = matrix(rnorm(d*(K-1),0,1), ncol = K-1, nrow = d); gamma_true # d by K-1
    if(equal_probs) gamma_true[1,] = c(-2.8,-2.5,-2.2,-1.3,-0.5) else gamma_true[1,] = c(1.8,0.5,0,0,0) 
    #gamma_true[1,] = c(1.8,0.5,0,0,0) # data probs not null 
    #gamma_true[1,] = c(-2.8,-2.5,-2.2,-1.3,-0.5) # equal probs not null
  }else{
    beta_true = matrix(0, ncol = K-1, nrow = q); beta_true # q by K-1
    gamma_true = matrix(0, ncol = K-1, nrow = d); gamma_true # d by K-1
    if(equal_probs) gamma_true[1,] = c(-1.8,-1.5,-1,-.5,-.3) else gamma_true[1,] = c(1.2,0.8,0.6,0.2,0.1)
    #gamma_true[1,] = c(1.2,0.8,0.6,0.2,0.1) # data probs null 
    #gamma_true[1,] = c(-1.8,-1.5,-1,-.5,-.3) # equal probs null
  }

  xbetaK_true = crossprod(t(x),beta_true); head(xbetaK_true) # n by K-1
  wgammaK_true = crossprod(t(w),gamma_true); head(wgammaK_true) # n by K-1
  psi_true = xbetaK_true + wgammaK_true

  # make pi
  pitildek = exp(psi_true)/(1+exp(psi_true));colMeans(pitildek)
  piik = matrix(0, n, K)
  piik[,1] = pitildek[,1]
  for(k in 2:(K-1)){
    piik[,k] = pitildek[,k]*(1-rowSums(matrix(piik[,(1:(k-1))],nrow=n,ncol = k-1)))
  }
  piik[,K] = 1-rowSums(piik)
  colMeans(piik)

  # simulate multinomial data
  # make sure each category has at least one known observation 
  repeat{
    y=t(sapply(1:n, FUN = function(i){
      return(rmultinom(n=1,size=1,prob=piik[i,]))
    }))
    if( !(0 %in% colMeans(y, na.rm = TRUE)) ) break
  }
  
  # which category 
  ymax = apply(y, 1, FUN = function(i){
    which(i==1)
  })
  
  allprobs = table(ymax) # save this for each simulation 
  #allprobs

  # save y
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
  
  # key to which ones are missing 
  ymiss = matrix(as.numeric(is.na(y)), nrow = n) # save ymiss 
  #length(which(is.na(y)))
  
  ycomplete_k = sapply(1:n, FUN = function(i){
    which(ycomplete[i,]==1)
  })
  
  y_save = y # save the missing pattern 
  
  #################################
  ### investigating missingness ###
  #################################
  
  # level of missing within each category 
  head(ymiss)
  colMeans(ymiss) 
  
  # # columns are categories 
  # ink = colSums(ycomplete) # how many truly assigned 
  # missk = colSums(ymiss) # how many missing
  # 
  # # of the missing data, how many are 1's 
  # num1 = c(sum(ymiss[,1][which(ycomplete_k==1)]),
  #                 sum(ymiss[,2][which(ycomplete_k==2)]),
  #                 sum(ymiss[,3][which(ycomplete_k==3)]),
  #                 sum(ymiss[,4][which(ycomplete_k==4)]),
  #                 sum(ymiss[,5][which(ycomplete_k==5)]),
  #                 sum(ymiss[,6][which(ycomplete_k==6)]))
  # 
  # num0 = c(sum(ymiss[,1][which(ycomplete_k!=1)]),
  #          sum(ymiss[,2][which(ycomplete_k!=2)]),
  #          sum(ymiss[,3][which(ycomplete_k!=3)]),
  #          sum(ymiss[,4][which(ycomplete_k!=4)]),
  #          sum(ymiss[,5][which(ycomplete_k!=5)]),
  #          sum(ymiss[,6][which(ycomplete_k!=6)]))
  # 
  # invest = rbind(ink, missk, num1, num0, t(num_meas_k))
  # 
  # 
  # df_invest = data.frame(cbind(invest, rowSums(invest))); df_invest
  # 
  # rownames(df_invest) = c("truly equal 1", "missing", "missing, equal 1", "missing, equal 0",
  #                         "true pos.", "true neg.", "false pos.", "false neg.")
  # 
  # colnames(df_invest) = c("cat 1", "cat 2", "cat 3", "cat 4", "cat 5", "cat 6", "total")
  # 
  # xtable(df_invest, digits = 0 )
  # 
  # 
  # save_y = y 
  
  ####################
  ### delete above ###
  ####################
  
  #################
  ### Fit Model ###
  #################
  
  fit = pgmultinom(niter = niter, priors = priors, y = y, ycomplete = ycomplete, x = x, w = w,
                   intercept = TRUE, beta_true = beta_true, gamma_true = gamma_true)
  
  
  if(miss_prob > 0){
    # complete case fit 
    y_cc = y[-whichmiss,]
    x_cc = x[-whichmiss,]
    w_cc = w[-whichmiss,]

    fit_cc = pgmultinom(niter = niter, priors = priors, y = y_cc, ycomplete = NULL, x = x_cc, w = w_cc,
                        intercept = TRUE, beta_true = beta_true, gamma_true = gamma_true)
  }

  
  
  ###################
  ### Imputations ###
  ###################
  
  if(miss_prob > 0){
    pscore = matrix(unlist(lapply((nburn + 1):niter, FUN = function(s){
      apply(fit$y.save[[s]][whichmiss,], 1, FUN = function(i){
        which(i==1)
      })
    })), nrow = niter-nburn, ncol = nmiss, byrow = TRUE)
    
    # each row is an iteration
    # each column is a case with missing outcome data 
    # values are the category assigned 
    
    # posterior probability of assignment to each cluster, propensity score 
    pprobs = t(sapply(1:nmiss, FUN = function(s){
      cat = pscore[,s]
      pcat = numeric(K)
      for(k in 1:K){
        pcat[k] = length(which(cat==k))/(niter-nburn)
      }
      return(pcat)
    }))
    
    # 0-1 loss
    # consider only the outcomes that are uncertain 
    # loss of true values will be 0 so divide by length(which(is.na(y)))
    #zero_one_loss = sum(abs(pprobs - ycomplete[whichmiss, ]))/length(which(is.na(y))) 
    #zero_one_loss = sum(abs(pprobs - ycomplete[whichmiss, ]))/length(y) 
    
    # 0-1 loss for each category 
    zero_one_loss = numeric(K)
    for(k in 1:K){
      zero_one_loss[k] = sum(abs(pprobs[,k] - ycomplete[whichmiss, k]))/length(which(ymiss[,k]==1))
    }

    # posterior mean classification evaluations 
    accuracy_mean = apply(fit$accuracy[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
    precision_mean = apply(fit$precision[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
    recall_mean = apply(fit$recall[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
    specificity_mean = apply(fit$specificity[(nburn+1):niter, ], 2, FUN = function(m) mean(m, na.rm = TRUE))
  
  }else{
    zero_one_loss = rep(NA, K)
    accuracy_mean = rep(NA, K)
    precision_mean = rep(NA, K)
    recall_mean = rep(NA, K)
    specificity_mean = rep(NA, K)
  }
  
  ########################################
  ### Regression coefficient estimates ###
  ########################################
  
  rmse_beta = mean(fit$rmse_beta[(nburn+1):niter])
  rmse_gamma = mean(fit$rmse_gamma[(nburn+1):niter])
  
  bias_beta = mean(fit$bias_beta[(nburn+1):niter])
  bias_gamma = mean(fit$bias_gamma[(nburn+1):niter]) 
  
  
  beta_quants = t(apply(fit$beta.vec[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), quantile(b, 0.975)))
  }))
  wid_beta = mean(beta_quants[,2]-beta_quants[,1]) # width of interval 
  cov_beta = mean(as.vector(beta_true[,-K]) > beta_quants[,1] & as.vector(beta_true[,-K]) < beta_quants[,2])
  
  gamma_quants = t(apply(fit$gamma.vec[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), quantile(b, 0.975)))
  }))
  wid_gamma = mean(gamma_quants[,2]-gamma_quants[,1])
  cov_gamma = mean(as.vector(gamma_true[,-d]) > gamma_quants[,1] & as.vector(gamma_true[,-d]) < gamma_quants[,2])

  
  # par(mfrow = c(3,5))
  # for(b in 1:15){
  #   plot(fit$beta.vec[,b], type = "l")
  #   abline(h = as.vector(beta_true)[b], col = "red", lwd = 3)
  # }

  #############################
  ### save table of results ###
  #############################
  
  df1 = data.frame(rbind(c(zero_one_loss, 
                          rmse_beta, bias_beta, wid_beta, cov_beta, 
                          rmse_gamma, bias_gamma, wid_gamma, cov_gamma,
                          accuracy_mean, precision_mean,
                          recall_mean, specificity_mean)))

  colnames(df1) = c("binary_loss1", "binary_loss2", "binary_loss3", "binary_loss4", "binary_loss5", "binary_loss6", 
             "RMSE_beta", "bias_beta", "ci_width_beta", "cov_beta", 
             "RMSE_gamma", "bias_gamma", "ci_width_gamma", "cov_gamma", 
             "accuracy1","accuracy2", "accuracy3","accuracy4","accuracy5","accuracy6",
             "precision1", "precision2", "precision3", "precision4", "precision5", "precision6", 
             "recall1","recall2","recall3","recall4","recall5","recall6", 
             "specificity1","specificity2","specificity3","specificity4","specificity5","specificity6")
  
  df_all = rbind(df_all, df1)
  
  
  ######################################################
  ### Complete Case regression coefficient estimates ###
  ######################################################
  

  if(miss_prob > 0){
    cc_rmse_beta = mean(fit_cc$rmse_beta[(nburn+1):niter])
    cc_rmse_gamma = mean(fit_cc$rmse_gamma[(nburn+1):niter])
    
    cc_bias_beta = mean(fit_cc$bias_beta[(nburn+1):niter])
    cc_bias_gamma = mean(fit_cc$bias_gamma[(nburn+1):niter]) 
    
    
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








