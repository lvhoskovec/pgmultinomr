
if ( requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("pgmultinomr")
}

library(Rcpp)
library(RcppArmadillo)
library(pgdraw)
library(mvnfast)

compileAttributes()
devtools::build()
devtools::install()