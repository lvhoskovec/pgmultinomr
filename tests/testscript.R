### script to test code in vignette ###
#######################################
rm(list = ls())
gc()

#################
### libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(pgdraw)
library(mvnfast)
compileAttributes()
devtools::build()
devtools::install()
library(pgmultinomr)
library(ggplot2)
library(cowplot)

#######################
### R markdown code ###
#######################

covX = matrix(c(1, -0.13, 0.02, -0.13, 1, -0.86, 0.02, -0.86, 1), 3, 3)

n = 1000
K = 6
p = 3
q = 5
dat = sim_dat(n=n, miss_prob = 0.5, K=K, p=p, q=q,
              covX = covX, allmiss = FALSE, null_scenario = FALSE, equal_probs = FALSE)

niter = 1000
fit = pgmultinom(niter=niter, priors=NULL, y = dat$y, ycomplete = dat$ycomplete, x = dat$x, w = dat$w, intercept = TRUE,
                 beta_true = dat$beta_true, gamma_true = dat$gamma_true)

nburn = niter/2
beta_quantiles = t(apply(fit$beta.vec[(nburn+1):niter,], 2, FUN = function(b){
  return(c(quantile(b, 0.025), mean(b), quantile(b, 0.975)))}))
namesX = c("exposure1", "exposure2", "exposure3")
namesY = c("cat1", "cat2", "cat3", "cat4", "cat5")
mix_df = cbind(data.frame(exp(beta_quantiles)), rep(namesX,K-1))
colnames(mix_df) = c("lwr", "mean", "upr", "exposure")
plot_list = list()

for(p in 1:length(namesX)){
  data_plot = mix_df[which(mix_df$exposure==namesX[p]),]
  limits = c(min(data_plot$lwr), max(data_plot$upr))
  # make a single exposure plot 
  mix_plot = ggplot() + geom_point(data = data_plot, mapping = aes(x = 1:(K-1), y = mean), size = 3) + 
    geom_errorbar(data = data_plot, mapping = aes(x = 1:(K-1), ymin = lwr, ymax = upr),width = 0.4) + 
    geom_hline(yintercept = 1, colour = "red") + 
    scale_x_discrete(name = "", breaks = factor(1:(K-1)), limits = factor(1:(K-1)), label = namesY) + 
    scale_y_continuous(name = "exp(beta)", limits = limits) + 
    theme(text = element_text(size = 15)) + ggtitle(namesX[p])
  # save plots in a list
  plot_list[[p]] = mix_plot
}

# plot all exposures 
plot_grid(
  plot_list[[1]] + theme(legend.position="none"),
  plot_list[[2]] + theme(legend.position="none"),
  plot_list[[3]] + theme(legend.position="none"),
  align = 'h',
  labels = "",
  hjust = -1,
  nrow = 1
)


len = 100
x.star = matrix(0, nrow = len, ncol = ncol(dat$x))
x.star[,1] = seq(-1,1,length.out = len) # sequence of exposure 1 values
x.star[,2] = quantile(dat$x[,2], 0.50) # set exposure 2 to a fixed percentile
x.star[,3] = quantile(dat$x[,3], 0.50) # set exposure 3 to a fixed percentile

# calculate the interactions
#x.star[,4] = x.star[,1]*x.star[,2]
#x.star[,5] = x.star[,1]*x.star[,3]
#x.star[,6] = x.star[,2]*x.star[,3]

# make x.base 
x.base = matrix(0, nrow = 1, ncol = ncol(dat$x))
x.base[,1] = mean(dat$x[,1]) # set exposure 1 to 0 
x.base[,2] = quantile(dat$x[,2], 0.50) # set exposure 2 to a fixed percentile
x.base[,3] = quantile(dat$x[,3], 0.50) # set exposure 3 to a fixed percentile

# calculate the interactions 
#x.base[,4] = x.base[,1]*x.base[,2]
#x.base[,5] = x.base[,1]*x.base[,3]
#x.base[,6] = x.base[,2]*x.base[,3]

# calculate OR(x.star, x.base)
odds_ratio = get_odds_ratio(x.star = x.star, x.base = x.base, 
                            beta_posterior = fit$beta, 
                            niter = 1000, nburn = 500, K = 6, 
                            refK = 3)
# make data frame of odds ratios
or_df = data.frame(xvals = x.star[,1], matrix(unlist(odds_ratio), nrow = len))
# set column names
or_names = "xvals"
for(k in 1:(K-1)){
  k_names = c("lwr", "mean", "upr")
  or_names = c(or_names,k_names)
}
colnames(or_df) = or_names
# make K-1 plots of odds ratio for each category
plot_or = list()
for(k in 1:(K-1)){
  # get data for category k 
  col_nums = 1 + (3*k-2):(3*k)
  df_k = or_df[,c(1,col_nums)]
  # plot OR for category k 
  plotk = ggplot(data = df_k) + ggtitle(paste0("cat",k)) + 
    geom_line(aes(x = xvals, y = mean), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr, ymax = upr), fill = "black",
                alpha = .4, show.legend = FALSE) +
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = "exposure1") + 
    scale_y_continuous(name = "estimated OR relative to mean") +
    theme(text = element_text(size = 15),
          axis.title.y = element_text(size = 12))
  plot_or[[k]] = plotk
}

plot_grid(
  plot_or[[1]] + theme(legend.position="none"),
  plot_or[[2]] + theme(legend.position="none"),
  plot_or[[3]] + theme(legend.position="none"),
  plot_or[[4]] + theme(legend.position="none"),
  plot_or[[5]] + theme(legend.position="none"),
  align = 'h',
  labels = "",
  hjust = -1,
  nrow = 2
)



