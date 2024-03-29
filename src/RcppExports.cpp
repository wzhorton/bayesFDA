// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bs_even
arma::mat bs_even(arma::vec time, int nk);
RcppExport SEXP _bayesFDA_bs_even(SEXP timeSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(bs_even(time, nk));
    return rcpp_result_gen;
END_RCPP
}
// g2g_fda
List g2g_fda(arma::mat curves, arma::vec lasts, arma::vec time, int p, int niter, int nburn);
RcppExport SEXP _bayesFDA_g2g_fda(SEXP curvesSEXP, SEXP lastsSEXP, SEXP timeSEXP, SEXP pSEXP, SEXP niterSEXP, SEXP nburnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type curves(curvesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lasts(lastsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    rcpp_result_gen = Rcpp::wrap(g2g_fda(curves, lasts, time, p, niter, nburn));
    return rcpp_result_gen;
END_RCPP
}
// cov_fda
List cov_fda(arma::mat curves, arma::mat X, arma::vec time, int p, int niter, int nburn);
RcppExport SEXP _bayesFDA_cov_fda(SEXP curvesSEXP, SEXP XSEXP, SEXP timeSEXP, SEXP pSEXP, SEXP niterSEXP, SEXP nburnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type curves(curvesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_fda(curves, X, time, p, niter, nburn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesFDA_bs_even", (DL_FUNC) &_bayesFDA_bs_even, 2},
    {"_bayesFDA_g2g_fda", (DL_FUNC) &_bayesFDA_g2g_fda, 6},
    {"_bayesFDA_cov_fda", (DL_FUNC) &_bayesFDA_cov_fda, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesFDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
