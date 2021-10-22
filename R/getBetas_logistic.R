#' Get summary of exposure regression coefficient estimates in logistic model 
#'
#' @param fit1 object of type 'pglogistic'
#' @param nburn number of burn-in iterations
#' @param niter number of total iterations
#'
#' @return summary of posterior distribution of regression coefficients
#' @export
#'

getBetas_logistic = function(fit1, nburn, niter){
  
  # beta estimates, posterior mean, 95% credible interval 
  beta_quants = t(apply(fit1$beta[(nburn+1):niter,], 2, FUN = function(b){
    return(c(quantile(b, 0.025), mean(b), quantile(b, 0.975)))
  }))
  sig_beta = apply(beta_quants, 1, FUN = function(x){
    return(x[3] < 0 | x[1] > 0)
  })
  beta_mean = beta_quants[,2]
  
  summary_list = list(beta_quants = beta_quants, which_sig_beta = which(sig_beta), 
                      exp_beta_quants = exp(beta_quants))
  

}
