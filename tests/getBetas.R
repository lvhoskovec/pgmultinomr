#' Get summary of exposure regression coefficient estimates
#'
#' @param fit1 object of type 'pgmultinom'
#' @param nburn number of burn-in iterations
#' @param niter total number of iterations 
#' @param dimX dimension of exposure data 
#' @param K dimension of outcome data 
#' @param namesX names of exposure data in order of columns of x
#' @param namesY names of outcome categories in order of y for k = 2,...,K
#'
#' @return summary of posterior distribution of regression coefficients
#' @export
#'

getBetas = function(fit1, nburn, niter, dimX, K, namesX = NULL, namesY = NULL){
  
  # beta estimates, posterior mean, 95% credible interval 
  beta_quants = t(apply(fit1$beta.vec[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), mean(b), quantile(b, 0.975)))
  }))
  beta_quants
  exp(beta_quants)
  sig_beta = apply(beta_quants, 1, FUN = function(x){
    return(x[3] < 0 | x[1] > 0)
  })
  which(sig_beta)
  beta_mean = beta_quants[,2]
  beta_mean[which(sig_beta)]
  
  # posterior mean
  beta_mean_df = data.frame(matrix(beta_mean, nrow = dimX, ncol = K-1))
  rownames(beta_mean_df) = namesX
  colnames(beta_mean_df) = namesY
  beta_mean_df
  exp(beta_mean_df)
  
  
  summary_list = list(posterior_mean_df = beta_mean_df, exp_post_mean_df = exp(beta_mean_df),
                      sig_beta = beta_mean[which(sig_beta)], exp_sig_beta = exp(beta_mean[which(sig_beta)]),
                      beta_quantiles = beta_quants, exp_beta_quantiles = exp(beta_quants),
                      which_sig_beta = which(sig_beta))
  
  
  return(summary_list)
  
  
  
}





