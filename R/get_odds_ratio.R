#' Get posterior distribution of odds ratio
#'
#' @param x.star x.star matrix for odds ratio
#' @param x.base x.base matrix for odds ratio
#' @param beta_posterior (list) posterior distribution of betas from model fit 
#' @param niter total number of MCMC iterations for inference 
#' @param nburn number of burn-in iterations
#' @param K number of outcome categories
#' @param refK index of reference outcome category
#'
#' @return posterior distribution of odds ratio 
#' @export 
#'

get_odds_ratio = function(x.star, x.base, beta_posterior, niter, nburn, K, refK){
  
  iter_seq = (nburn+1):niter
  n.star = nrow(x.star)
  # make space 
  piik_dist = list()
  for(k in 1:K){
    piik_dist[[k]] = matrix(NA, nrow = n.star, ncol = 1)
  }
  # P(Y = k|x = 0) incorporating distn of betas, an element of the list for each k 
  for(s in iter_seq){
    beta_iter = beta_posterior[[s]] # list 
    xbetaK_iter = crossprod(t(x.star),beta_iter); xbetaK_iter # n by K
    psi_iter = xbetaK_iter
    # make pi
    pitildek_iter = exp(psi_iter)/(1+exp(psi_iter))
    piik_iter = matrix(0, n.star, K-1)
    piik_iter[,1] = pitildek_iter[,1]
    for(k in 2:(K-1)){
      piik_iter[,k] = pitildek_iter[,k]*(1-rowSums(matrix(piik_iter[,(1:(k-1))],nrow=n.star,ncol = k-1)))
    }
    piik_iter = cbind(piik_iter, 1 - rowSums(piik_iter))
    for(k in 1:K){
      piik_dist[[k]] = cbind(piik_dist[[k]], piik_iter[,k])
    }
  }
  # remove the first NA column 
  for(k in 1:K){
    piik_dist[[k]] = piik_dist[[k]][,-1]
  }
  
  # odds relative to reference category
  odd_star = list()
  k_step = 1
  for(k in ((1:K)[-refK])){
    odd_star[[k_step]] = piik_dist[[k]]/piik_dist[[refK]]
    k_step=k_step+1
  }
  
  ### calculate distn of odds for x.base ###
  n.base = 1
  piik_dist = list()
  for(k in 1:K){
    piik_dist[[k]] = matrix(NA, nrow = n.base, ncol = 1)
  }
  for(s in iter_seq){
    beta_iter = beta_posterior[[s]] # list 
    xbetaK_iter = crossprod(t(x.base),beta_iter); xbetaK_iter # n by K
    psi_iter = xbetaK_iter
    # make pi
    pitildek_iter = exp(psi_iter)/(1+exp(psi_iter))
    piik_iter = matrix(0, n.base, K-1)
    piik_iter[,1] = pitildek_iter[,1]
    for(k in 2:(K-1)){
      piik_iter[,k] = pitildek_iter[,k]*(1-rowSums(matrix(piik_iter[,(1:(k-1))],nrow=n.base,ncol = k-1)))
    }
    piik_iter = cbind(piik_iter, 1 - rowSums(piik_iter))
    for(k in 1:K){
      piik_dist[[k]] = cbind(piik_dist[[k]], piik_iter[,k])
    }
  }
  # remove first column of NA 
  for(k in 1:K){
    piik_dist[[k]] = piik_dist[[k]][,-1]
  }
  
  odd_base = list()
  k_step = 1
  for(k in ((1:K)[-refK])){
    odd_base[[k_step]] = piik_dist[[k]]/piik_dist[[refK]]
    k_step=k_step+1
  }
  
  ## get posterior distn of odds ratio ##
  # P(Y = k|X = x)/P(Y = 0|X = x) # odd_star 
  # P(Y = k|X = 0)/P(Y = 0|X = 0) # odd_base 
  OR_summary = list()
  for(k in 1:(K-1)){
    OR_dist = matrix(unlist(lapply(1:n.star, FUN = function(i){
      odd_star[[k]][i,]/odd_base[[k]]
    })), ncol = (niter-nburn), nrow = n.star, byrow = TRUE)
    OR_summary[[k]] = t(apply(OR_dist, 1, FUN = function(x){
      c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
    }))
  }
  return(OR_summary)
}
