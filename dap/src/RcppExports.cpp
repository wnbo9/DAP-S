// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_log10_prior
double compute_log10_prior(IntegerVector mcfg, NumericVector pi_vec);
RcppExport SEXP _dap_compute_log10_prior(SEXP mcfgSEXP, SEXP pi_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type mcfg(mcfgSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_prior(mcfg, pi_vec));
    return rcpp_result_gen;
END_RCPP
}
// compute_log10_posterior
NumericVector compute_log10_posterior(NumericMatrix X, NumericVector y, IntegerMatrix cmfg_matrix, NumericVector pi_vec, NumericVector phi2_vec);
RcppExport SEXP _dap_compute_log10_posterior(SEXP XSEXP, SEXP ySEXP, SEXP cmfg_matrixSEXP, SEXP pi_vecSEXP, SEXP phi2_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type cmfg_matrix(cmfg_matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi2_vec(phi2_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_posterior(X, y, cmfg_matrix, pi_vec, phi2_vec));
    return rcpp_result_gen;
END_RCPP
}
// pir
NumericMatrix pir(NumericMatrix mat, double threshold);
RcppExport SEXP _dap_pir(SEXP matSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(pir(mat, threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_compute_log10_prior", (DL_FUNC) &_dap_compute_log10_prior, 2},
    {"_dap_compute_log10_posterior", (DL_FUNC) &_dap_compute_log10_posterior, 5},
    {"_dap_pir", (DL_FUNC) &_dap_pir, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
