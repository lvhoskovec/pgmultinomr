#####################################################
### functions for analyzing data analysis results ###
#####################################################


## get the posterior distribution of the odds ratio ##
getOR_dist = function(x.star, x.base, fit, niter, nburn, incr = FALSE){
  
  iter_seq = (nburn+1):niter
  ## calculate distn of odds for x.star ##
  n.star = nrow(x.star)
  # make space 
  piik_dist = list()
  for(k in 1:K){
    piik_dist[[k]] = matrix(NA, nrow = n.star, ncol = 1)
  }
  # P(Y = k|x = 0) incorporating distn of betas, an element of the list for each k 
  for(s in iter_seq){
    beta_iter = fit$beta[[s]] # list 
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
  odd_star = list()
  if(incr){
    # icnremental OR
    odd_star[[1]] = piik_dist[[1]]/piik_dist[[2]]   # sympt wrt asympt 
    odd_star[[2]] = piik_dist[[3]]/(piik_dist[[2]]+piik_dist[[1]]) # hosp wrt sympt and asympt
    odd_star[[3]] = piik_dist[[4]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]])
    odd_star[[4]] = piik_dist[[5]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]]+piik_dist[[4]])
    odd_star[[5]] = piik_dist[[6]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]]+piik_dist[[4]]+piik_dist[[5]])
  }else{
    # divide each category by cat 2 to get odds relative to asymptomatic 
    kk = 1
    for(k in c(1,3,4,5,6)){
      odd_star[[kk]] = piik_dist[[k]]/piik_dist[[2]]
      kk=kk+1
    }
  }
  
  ## calculate distn of odds for x.base ##
  n.base = nrow(x.base)
  # make space
  piik_dist = list()
  for(k in 1:K){
    piik_dist[[k]] = matrix(NA, nrow = n.base, ncol = 1)
  }
  # P(Y = k|x = 0) incorporating distn of betas, an element of the list for each k 
  for(s in iter_seq){
    beta_iter = fit$beta[[s]] # list 
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
  if(incr){
    # icnremental OR
    odd_base[[1]] = piik_dist[[1]]/piik_dist[[2]]   # sympt wrt asympt 
    odd_base[[2]] = piik_dist[[3]]/(piik_dist[[2]]+piik_dist[[1]]) # hosp wrt sympt and asympt
    odd_base[[3]] = piik_dist[[4]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]])
    odd_base[[4]] = piik_dist[[5]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]]+piik_dist[[4]])
    odd_base[[5]] = piik_dist[[6]]/(piik_dist[[2]]+piik_dist[[1]]+piik_dist[[3]]+piik_dist[[4]]+piik_dist[[5]])
  }else{
    # divide each category by cat 2 to get odds relative to asymptomatic 
    kk = 1
    for(k in c(1,3,4,5,6)){
      odd_base[[kk]] = piik_dist[[k]]/piik_dist[[2]]
      kk=kk+1
    }
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

## make x matrix of values to calculate odds ratio for mixture with 3 pollutants and interactions ##
make_xstar = function(len=100, x_seq, x, quant = 0.25, move_exp = 1){
  
  x.star = matrix(0, nrow = len, ncol = 6)
  # changes in this exposure 
  x.star[,move_exp] = x_seq 
  # holding these exposures at quant 
  others = (1:3)[-move_exp]
  
  if(length(quant) == 1){
    x.star[,others] = matrix(apply(x[,others], 2, FUN = function(xx) quantile(xx, quant)), 
                             nrow = len, ncol = 2, byrow = TRUE)
  }else if(length(quant)==2){
    others1 = others[1]
    others2 = others[2]
    x.star[,others1] = quantile(x[,others1], quant[1])
    x.star[,others2] = quantile(x[,others2], quant[2])
  }


  # calculate the interactions 
  x.star[,4] = x.star[,1]*x.star[,2]
  x.star[,5] = x.star[,1]*x.star[,3]
  x.star[,6] = x.star[,2]*x.star[,3]
  return(x.star)
}



## make the odds ratio plots ##
make_oddsRatio_interval_plots = function(len, x_seq = seq(-2,2,length.out = 100),
                                         x, expdata, quant, move_exp = 1,
                                         fit, niter, nburn,
                                         exp_name = "PM", incr = FALSE,
                                         model = "mixture", 
                                         outcome_labs = c("symptomatic",
                                                          "hospitalized",
                                                          "ICU",
                                                          "ventilator",
                                                          "death"),
                                         ylab = "Odds Ratio"){
  
  
  
  if(model == "mixture"){
    x.star = make_xstar(len = len, x_seq = x_seq,
                        x = x, quant = quant, move_exp = move_exp)
    x.base = make_xstar(len = 1, x_seq = 0,
                        x = x, quant = quant, move_exp = move_exp)
    odds_dist = getOR_dist(x.star = x.star, x.base = x.base, 
                           fit = fit, niter = niter, nburn = nburn, incr = incr)
  }else if(model == "single"){
    x.star = matrix(x_seq, ncol = 1)
    x.base = matrix(0, nrow = 1, ncol = 1)
    
    odds_dist = getOR_dist(x.star = x.star, x.base = x.base, 
               fit = fit, niter = niter, nburn = nburn, incr = incr)
  }

  
  
  odds_interval_df = data.frame(xvals = x.star[,move_exp],
                                lwr1 = odds_dist[[1]][,1],
                                mean1 = odds_dist[[1]][,2],
                                upr1 = odds_dist[[1]][,3],
                                lwr2 = odds_dist[[2]][,1],
                                mean2 = odds_dist[[2]][,2],
                                upr2 = odds_dist[[2]][,3],
                                lwr3 = odds_dist[[3]][,1],
                                mean3 = odds_dist[[3]][,2],
                                upr3 = odds_dist[[3]][,3],
                                lwr4 = odds_dist[[4]][,1],
                                mean4 = odds_dist[[4]][,2],
                                upr4 = odds_dist[[4]][,3],
                                lwr5 = odds_dist[[5]][,1],
                                mean5 = odds_dist[[5]][,2],
                                upr5 = odds_dist[[5]][,3])
  
  
  
  # x labels

  xbreaks = round(seq(min(x_seq), max(x_seq), length.out = 9),1)
  iqr = IQR(expdata[,move_exp])
  mn = mean(expdata[,move_exp])
  xlabs = round(xbreaks*iqr + mn,1)
  
  # combine plot mean only
  mean_max = max(odds_interval_df$mean1, odds_interval_df$mean2, odds_interval_df$mean3, 
                 odds_interval_df$mean4, odds_interval_df$mean5)
  mean_plot = ggplot(data = odds_interval_df) + geom_line(aes(x = xvals, y = mean1, col = factor(1)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean2, col = factor(2)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean3, col = factor(3)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean4, col = factor(4)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean5, col = factor(5)), size = 1.5) + 
    labs(col = "Category") + 
    scale_color_discrete(labels = outcome_labs) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,mean_max), breaks = round(seq(0,mean_max,length.out = 10),1)) 
  
  mean_plot
  
  # combine plot with intervals
  upr_max = max(odds_interval_df$upr1, odds_interval_df$upr2, odds_interval_df$upr3, 
                odds_interval_df$upr4, odds_interval_df$upr5)
  mean_ci_plot = ggplot(data = odds_interval_df) + 
    geom_line(aes(x = xvals, y = mean1, col = factor(1)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean2, col = factor(2)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean3, col = factor(3)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean4, col = factor(4)), size = 1.5) + 
    geom_line(aes(x = xvals, y = mean5, col = factor(5)), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr1, ymax = upr1), fill = "red", col = "red", size = 2, 
                alpha = .4, show.legend = FALSE) +
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr2, ymax = upr2), fill = "yellow",col = "yellow", size = 2, 
                alpha = .4, show.legend = FALSE) +
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr3, ymax = upr3), fill = "green",col = "green", size = 2, 
                alpha = .4, show.legend = FALSE) +
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr4, ymax = upr4), fill = "blue",col = "blue", size = 2, 
                alpha = .4, show.legend = FALSE) +
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr5, ymax = upr5), fill = "pink",col = "pink", size = 2, 
                alpha = .4, show.legend = FALSE) +
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) +
    labs(col = "Category") + 
    scale_color_discrete(labels = outcome_labs) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,upr_max), breaks = round(seq(0,upr_max,by=0.1),1)) 
  
  # symptomatic
  s_max = max(odds_interval_df$upr1)
  if(s_max <= 1.2) {
    s_by = 0.1
  } else if(s_max <=2.6){
    s_by = 0.2
  } else s_by = 0.4
  s_breaks = round(seq(0,s_max,by=s_by),1)
  
  s_plot = ggplot(data = odds_interval_df) + ggtitle("symptomatic") + geom_line(aes(x = xvals, y = mean1), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr1, ymax = upr1), fill = "black",
                alpha = .4, show.legend = FALSE) +
    labs(col = "Category") + 
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,s_max), breaks = s_breaks) +
    theme(text = element_text(size = 20))
  
  # hospitalized 
  h_max = max(odds_interval_df$upr2)
  if(h_max <= 1.2) {
    h_by = 0.1
  } else if(h_max <=2.6){
    h_by = 0.2
  } else h_by = 0.4
  h_breaks = round(seq(0,h_max,by=h_by),1)
  
  h_plot = ggplot(data = odds_interval_df) + ggtitle("hospitalized") + geom_line(aes(x = xvals, y = mean2), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr2, ymax = upr2), fill = "black",
                alpha = .4, show.legend = FALSE) +
    labs(col = "Category") + 
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,h_max), breaks = h_breaks) +
    theme(text = element_text(size = 20))
  
  # icu
  i_max = max(odds_interval_df$upr3)
  if(i_max <= 1.2) {
    i_by = 0.1
  } else if(i_max <=2.6){
    i_by = 0.2
  } else i_by = 0.4
  i_breaks = round(seq(0,i_max,by=i_by),1)
  
  
  i_plot = ggplot(data = odds_interval_df) + ggtitle("ICU") + geom_line(aes(x = xvals, y = mean3), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr3, ymax = upr3), fill = "black",
                alpha = .4, show.legend = FALSE) +
    labs(col = "Category") + 
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,i_max), breaks = i_breaks) +
    theme(text = element_text(size = 20))
  
  # vent
  v_max = max(odds_interval_df$upr4)
  if(v_max <= 1.2) {
    v_by = 0.1
  } else if(v_max <=2.6){
    v_by = 0.2
  } else v_by = 0.4
  v_breaks = round(seq(0,v_max,by=v_by),1)
  
  v_plot = ggplot(data = odds_interval_df) + ggtitle("ventilator") + geom_line(aes(x = xvals, y = mean4), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr4, ymax = upr4), fill = "black",
                alpha = .4, show.legend = FALSE) +
    labs(col = "Category") + 
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,v_max), breaks = v_breaks) +
    theme(text = element_text(size = 20))
  
  # death 
  d_max = max(odds_interval_df$upr5)
  if(d_max <= 1.2) {
    d_by = 0.1
  } else if(d_max <=2.6){
    d_by = 0.2
  } else d_by = 0.4
  d_breaks = round(seq(0,d_max,by=d_by),1)
  
  d_plot = ggplot(data = odds_interval_df) + ggtitle("death") + geom_line(aes(x = xvals, y = mean5), size = 1.5) + 
    geom_ribbon(mapping = 
                  aes(x = xvals, ymin = lwr5, ymax = upr5), fill = "black",
                alpha = .4, show.legend = FALSE) +
    labs(col = "Category") + 
    geom_line(aes(x = xvals, y = 1), linetype = "dotted", size = 1.5) + 
    scale_x_continuous(name = exp_name, breaks = xbreaks, labels = xlabs) + 
    scale_y_continuous(name = ylab, limits = c(0,d_max), breaks = d_breaks) +
    theme(text = element_text(size = 20))
  
  plot_list = list(mean_plot, 
                   mean_ci_plot, 
                   s_plot, 
                   h_plot, 
                   i_plot,
                   v_plot,
                   d_plot)
  
  return(plot_list)
  
  
}

# plot regression coefficients 
plotBetas = function(beta_quantiles, K, namesX, namesY, ylab = "exp(beta)"){
  
  
  exp_names = namesX
  mix_df = cbind(data.frame(exp(beta_quantiles)), rep(exp_names,K-1))
  colnames(mix_df) = c("lwr", "mean", "upr", "exposure")
  
  plot_list = list()
  
  for(p in 1:length(namesX)){
    
    plot_name = paste0("plot",p)
    data_plot = mix_df[which(mix_df$exposure==namesX[p]),]
    
    rg = seq(min(data_plot$lwr)-0.1,max(data_plot$upr)+0.1, by = 0.1)
    if(length(rg)>10) sp = .2 else sp = .1
    
    breaks = round(seq(min(data_plot$lwr)-0.1,max(data_plot$upr)+0.1, by = sp),1)
    limits = c(min(data_plot$lwr), max(data_plot$upr))
    
    
    mix_plot = ggplot() + geom_point(data = data_plot, mapping = aes(x = 1:(K-1), y = mean), size = 3) + 
      geom_errorbar(data = data_plot, mapping = aes(x = 1:(K-1), ymin = lwr, ymax = upr),width = 0.4) + 
      geom_hline(yintercept = 1, colour = "red") + 
      scale_x_discrete(name = "", breaks = factor(1:(K-1)), limits = factor(1:(K-1)), label = namesY) + 
      scale_y_continuous(name = ylab, breaks = breaks, limits = limits) + 
      theme(text = element_text(size = 15)) + ggtitle(namesX[p])
    
    plot_list[[p]] = mix_plot
    
  }
  
  return(plot_list)
  
}



