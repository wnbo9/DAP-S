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
List compute_log10_posterior(NumericMatrix X, NumericVector y, NumericMatrix cmfg_matrix, NumericVector pi_vec, NumericVector phi2_vec);
RcppExport SEXP _dap_compute_log10_posterior(SEXP XSEXP, SEXP ySEXP, SEXP cmfg_matrixSEXP, SEXP pi_vecSEXP, SEXP phi2_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cmfg_matrix(cmfg_matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi2_vec(phi2_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_posterior(X, y, cmfg_matrix, pi_vec, phi2_vec));
    return rcpp_result_gen;
END_RCPP
}
// dap_main
List dap_main(NumericMatrix X, NumericVector y, NumericMatrix matrix, double threshold, NumericVector prior_weights, NumericVector phi2_vec, double r2_threshold, double coverage);
RcppExport SEXP _dap_dap_main(SEXP XSEXP, SEXP ySEXP, SEXP matrixSEXP, SEXP thresholdSEXP, SEXP prior_weightsSEXP, SEXP phi2_vecSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_weights(prior_weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi2_vec(phi2_vecSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    rcpp_result_gen = Rcpp::wrap(dap_main(X, y, matrix, threshold, prior_weights, phi2_vec, r2_threshold, coverage));
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
// get_sc
List get_sc(const NumericMatrix& X, const NumericMatrix& mat, const NumericMatrix& cmfg_mat, const NumericVector& posterior_prob, const CharacterVector& col_names, double r2_threshold, double coverage);
RcppExport SEXP _dap_get_sc(SEXP XSEXP, SEXP matSEXP, SEXP cmfg_matSEXP, SEXP posterior_probSEXP, SEXP col_namesSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cmfg_mat(cmfg_matSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type posterior_prob(posterior_probSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type col_names(col_namesSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sc(X, mat, cmfg_mat, posterior_prob, col_names, r2_threshold, coverage));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_compute_log10_prior", (DL_FUNC) &_dap_compute_log10_prior, 2},
    {"_dap_compute_log10_posterior", (DL_FUNC) &_dap_compute_log10_posterior, 5},
    {"_dap_dap_main", (DL_FUNC) &_dap_dap_main, 8},
    {"_dap_pir", (DL_FUNC) &_dap_pir, 2},
    {"_dap_get_sc", (DL_FUNC) &_dap_get_sc, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
