#' Simulate data to fit pgmultinom 
#'
#' @param n sample size
#' @param miss_prob proportion of cases with missing outcome data 
#' @param K number of outcome categories 
#' @param p number of exposures
#' @param q number of covariates 
#' @param covX covariance of exposure data for simulation 
#' @param allmiss logical; if TRUE then all outcomes are missing for any case with missing outcome data, default is FALSE 
#' @param null_scenario logical; if TRUE then exposures and covariates have no an effect on outcome, default is FALSE
#' @param equal_probs logical; if TRUE the simulation setting is equal probabilities, if FALSE the simulation setting is data probabilities
#'
#' @importFrom mvnfast rmvn
#' @importFrom stats rnorm quantile rmultinom
#' @importFrom utils head
#'
#' @return summary of simulation results
#' @export
#'
#'
#'
sim_dat = function(n=1000, miss_prob = 0, K=6, p = 3, q = 5,
                     covX = diag(3), allmiss = FALSE, null_scenario = FALSE, equal_probs = FALSE){
  
  # simulate x data with covariance structure 
  p = dim(covX)[1]
  x = rmvn(n, mu = rep(0,p), sigma = covX) # change the mean to -1 
  
  # independent covariates
  q = 5
  w_cov = rmvn(n, mu = rep(0,q), sigma = diag(q))
  w = cbind(1, w_cov)
  q = ncol(w)
  
  # stick-breaking simulation #
  if(!null_scenario){
    beta_true = matrix(rnorm(p*(K-1),0,1), ncol = K-1, nrow = p) # p by K-1
    gamma_true = matrix(rnorm(q*(K-1),0,1), ncol = K-1, nrow = q) # q by K-1
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
  pitildek = exp(psi_true)/(1+exp(psi_true))
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
  
  # holdout set
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
  
  ldata = list(y = y, ycomplete = ycomplete, x = x, w = w, beta_true, gamma_true)
  return(ldata)
  
}
