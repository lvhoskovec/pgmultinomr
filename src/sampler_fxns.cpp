#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
List update_beta(int k, arma::vec ok, arma::vec yk, arma::mat xk, arma::vec wgam, 
                 List betaVarInverse, List betaMean){
  
  arma::vec zk = yk/ok - wgam; // matrix with 1 column 
  
  int ncol = xk.n_cols; 
  arma::mat ok_rep = ok;
  for(int j = 1; j < ncol; ++j){
    ok_rep = join_horiz(ok_rep,ok); // repeat ok ncol times for ncol columns of omegas
  }
  arma::mat xome = (xk%ok_rep).t(); // x times omega
  arma::mat Vk = inv(as<arma::mat>(betaVarInverse[k]) + xome*xk);
  arma::mat mk = Vk*(as<arma::mat>(betaVarInverse[k])*as<arma::colvec>(betaMean[k]) + xome*zk);
  List params(2);
  params[0]=mk;
  params[1]=Vk;
  return(params);
}

// [[Rcpp::export]]
List update_gamma(int k, arma::vec ok, arma::vec yk, arma::mat wk, arma::vec xbet, 
                 List gammaVarInverse, List gammaMean){
  
  arma::vec zk = yk/ok - xbet; // matrix with 1 column 
  
  int ncol = wk.n_cols; 
  arma::mat ok_rep = ok;
  for(int j = 1; j < ncol; ++j){
    ok_rep = join_horiz(ok_rep,ok); // repeat ok ncol times for ncol columns of omegas
  }
  arma::mat wome = (wk%ok_rep).t(); // x times omega
  arma::mat Vk = inv(as<arma::mat>(gammaVarInverse[k]) + wome*wk);
  arma::mat mk = Vk*(as<arma::mat>(gammaVarInverse[k])*as<arma::colvec>(gammaMean[k]) + wome*zk);
  List params(2);
  params[0]=mk;
  params[1]=Vk;
  return(params);
}


// [[Rcpp::export]]
List update_beta_logistic(arma::vec omega, arma::vec kappa, arma::mat x, arma::vec wgamma, 
                          arma::mat betaVarInverse, arma::vec betaMean){
  
  
  
  arma::vec z = kappa/omega - wgamma; // matrix with 1 column 
  
  int ncol = x.n_cols; 
  arma::mat o_rep = omega;
  for(int j = 1; j < ncol; ++j){
    o_rep = join_horiz(o_rep,omega); // repeat ok ncol times for ncol columns of omegas
  }
  arma::mat xome = (x%o_rep).t(); // x times omega
  arma::mat V = inv(betaVarInverse + xome*x);
  arma::mat m = V*(betaVarInverse*betaMean + xome*z);
  List params(2);
  params[0]=m;
  params[1]=V;
  return(params);
}

// [[Rcpp::export]]
List update_gamma_logistic(arma::vec omega, arma::vec kappa, arma::mat w, arma::vec xbeta, 
                  arma::mat gammaVarInverse, arma::vec gammaMean){
  
  arma::vec z = kappa/omega - xbeta; // matrix with 1 column 
  
  int ncol = w.n_cols; 
  arma::mat o_rep = omega;
  for(int j = 1; j < ncol; ++j){
    o_rep = join_horiz(o_rep,omega); // repeat ok ncol times for ncol columns of omegas
  }
  arma::mat wome = (w%o_rep).t(); // x times omega
  arma::mat V = inv(gammaVarInverse + wome*w);
  arma::mat m = V*(gammaVarInverse*gammaMean + wome*z);
  List params(2);
  params[0]=m;
  params[1]=V;
  return(params);
}




// [[Rcpp::export]]
arma::mat make_xbeta(arma::mat x, arma::mat beta){
  
  arma::mat xbetak = x * beta; 
  return(xbetak); 
  
}






