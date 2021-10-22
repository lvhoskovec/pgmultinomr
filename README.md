
# pgmultinomr

Fit Bayesian Categorical Regression Model with Polya-Gamma Augmentation

This folder contains code for the R package pgmultinomr.

The main purpose of this package is to fit a Bayesian categorical regression model using polya-gamma data augmentation for efficient Gibbs sampling. The model is based on the stick-breaking representation of the multinomial distribution. The method permits partially to fully missing outcome data and imputes values for the missing observations. We provide code for post-processing and visualizing the results.

This package relies on the following R packages. Install these packages in the R environment by using the install.packages("") command.

Rcpp

RcppArmadillo

mvnfast

pgdraw

Then you can install psbpHMM by running the following lines in the R console:

library(devtools)

install_github( "lvhoskovec/pgmultinomr", build_vignettes = TRUE)

library(pgmultinomr)

vignette("pgmultinomr_vignette")
