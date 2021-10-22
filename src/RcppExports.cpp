// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// update_omega
NumericMatrix update_omega(int n, int K, List whichik, NumericMatrix xbetak, NumericMatrix zgammak);
RcppExport SEXP _pgmultinomr_update_omega(SEXP nSEXP, SEXP KSEXP, SEXP whichikSEXP, SEXP xbetakSEXP, SEXP zgammakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< List >::type whichik(whichikSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xbetak(xbetakSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zgammak(zgammakSEXP);
    rcpp_result_gen = Rcpp::wrap(update_omega(n, K, whichik, xbetak, zgammak));
    return rcpp_result_gen;
END_RCPP
}
// fillMatrix
NumericMatrix fillMatrix(NumericVector a);
RcppExport SEXP _pgmultinomr_fillMatrix(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(fillMatrix(a));
    return rcpp_result_gen;
END_RCPP
}
// return_pgdraw
NumericVector return_pgdraw(NumericVector b, NumericVector c);
RcppExport SEXP _pgmultinomr_return_pgdraw(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(return_pgdraw(b, c));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pgdraw
NumericVector rcpp_pgdraw(NumericVector b, NumericVector c);
RcppExport SEXP _pgmultinomr_rcpp_pgdraw(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pgdraw(b, c));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _pgmultinomr_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _pgmultinomr_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _pgmultinomr_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _pgmultinomr_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// update_beta
List update_beta(int k, arma::vec ok, arma::vec yk, arma::mat xk, arma::vec wgam, List betaVarInverse, List betaMean);
RcppExport SEXP _pgmultinomr_update_beta(SEXP kSEXP, SEXP okSEXP, SEXP ykSEXP, SEXP xkSEXP, SEXP wgamSEXP, SEXP betaVarInverseSEXP, SEXP betaMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ok(okSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yk(ykSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xk(xkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wgam(wgamSEXP);
    Rcpp::traits::input_parameter< List >::type betaVarInverse(betaVarInverseSEXP);
    Rcpp::traits::input_parameter< List >::type betaMean(betaMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta(k, ok, yk, xk, wgam, betaVarInverse, betaMean));
    return rcpp_result_gen;
END_RCPP
}
// update_gamma
List update_gamma(int k, arma::vec ok, arma::vec yk, arma::mat wk, arma::vec xbet, List gammaVarInverse, List gammaMean);
RcppExport SEXP _pgmultinomr_update_gamma(SEXP kSEXP, SEXP okSEXP, SEXP ykSEXP, SEXP wkSEXP, SEXP xbetSEXP, SEXP gammaVarInverseSEXP, SEXP gammaMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ok(okSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yk(ykSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type wk(wkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xbet(xbetSEXP);
    Rcpp::traits::input_parameter< List >::type gammaVarInverse(gammaVarInverseSEXP);
    Rcpp::traits::input_parameter< List >::type gammaMean(gammaMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(update_gamma(k, ok, yk, wk, xbet, gammaVarInverse, gammaMean));
    return rcpp_result_gen;
END_RCPP
}
// update_beta_logistic
List update_beta_logistic(arma::vec omega, arma::vec kappa, arma::mat x, arma::vec wgamma, arma::mat betaVarInverse, arma::vec betaMean);
RcppExport SEXP _pgmultinomr_update_beta_logistic(SEXP omegaSEXP, SEXP kappaSEXP, SEXP xSEXP, SEXP wgammaSEXP, SEXP betaVarInverseSEXP, SEXP betaMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wgamma(wgammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type betaVarInverse(betaVarInverseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betaMean(betaMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta_logistic(omega, kappa, x, wgamma, betaVarInverse, betaMean));
    return rcpp_result_gen;
END_RCPP
}
// update_gamma_logistic
List update_gamma_logistic(arma::vec omega, arma::vec kappa, arma::mat w, arma::vec xbeta, arma::mat gammaVarInverse, arma::vec gammaMean);
RcppExport SEXP _pgmultinomr_update_gamma_logistic(SEXP omegaSEXP, SEXP kappaSEXP, SEXP wSEXP, SEXP xbetaSEXP, SEXP gammaVarInverseSEXP, SEXP gammaMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xbeta(xbetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gammaVarInverse(gammaVarInverseSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gammaMean(gammaMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(update_gamma_logistic(omega, kappa, w, xbeta, gammaVarInverse, gammaMean));
    return rcpp_result_gen;
END_RCPP
}
// make_xbeta
arma::mat make_xbeta(arma::mat x, arma::mat beta);
RcppExport SEXP _pgmultinomr_make_xbeta(SEXP xSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(make_xbeta(x, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pgmultinomr_update_omega", (DL_FUNC) &_pgmultinomr_update_omega, 5},
    {"_pgmultinomr_fillMatrix", (DL_FUNC) &_pgmultinomr_fillMatrix, 1},
    {"_pgmultinomr_return_pgdraw", (DL_FUNC) &_pgmultinomr_return_pgdraw, 2},
    {"_pgmultinomr_rcpp_pgdraw", (DL_FUNC) &_pgmultinomr_rcpp_pgdraw, 2},
    {"_pgmultinomr_rcpparma_hello_world", (DL_FUNC) &_pgmultinomr_rcpparma_hello_world, 0},
    {"_pgmultinomr_rcpparma_outerproduct", (DL_FUNC) &_pgmultinomr_rcpparma_outerproduct, 1},
    {"_pgmultinomr_rcpparma_innerproduct", (DL_FUNC) &_pgmultinomr_rcpparma_innerproduct, 1},
    {"_pgmultinomr_rcpparma_bothproducts", (DL_FUNC) &_pgmultinomr_rcpparma_bothproducts, 1},
    {"_pgmultinomr_update_beta", (DL_FUNC) &_pgmultinomr_update_beta, 7},
    {"_pgmultinomr_update_gamma", (DL_FUNC) &_pgmultinomr_update_gamma, 7},
    {"_pgmultinomr_update_beta_logistic", (DL_FUNC) &_pgmultinomr_update_beta_logistic, 6},
    {"_pgmultinomr_update_gamma_logistic", (DL_FUNC) &_pgmultinomr_update_gamma_logistic, 6},
    {"_pgmultinomr_make_xbeta", (DL_FUNC) &_pgmultinomr_make_xbeta, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pgmultinomr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}