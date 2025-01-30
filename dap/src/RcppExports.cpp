// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_log10_prior
double compute_log10_prior(const std::vector<int>& mcfg, NumericVector pi_vec);
RcppExport SEXP _dap_compute_log10_prior(SEXP mcfgSEXP, SEXP pi_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type mcfg(mcfgSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_prior(mcfg, pi_vec));
    return rcpp_result_gen;
END_RCPP
}
// compute_log10_posterior
List compute_log10_posterior(NumericMatrix X, NumericVector y, NumericMatrix cmfg_matrix, NumericMatrix single_matrix, NumericVector pi_vec, NumericMatrix phi2_mat);
RcppExport SEXP _dap_compute_log10_posterior(SEXP XSEXP, SEXP ySEXP, SEXP cmfg_matrixSEXP, SEXP single_matrixSEXP, SEXP pi_vecSEXP, SEXP phi2_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cmfg_matrix(cmfg_matrixSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type single_matrix(single_matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi2_mat(phi2_matSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_posterior(X, y, cmfg_matrix, single_matrix, pi_vec, phi2_mat));
    return rcpp_result_gen;
END_RCPP
}
// dap_main
List dap_main(NumericMatrix X, NumericVector y, NumericMatrix matrix, NumericVector prior_weights, double r2_threshold, double coverage, NumericMatrix phi2_mat, bool exclusive, double pir_threshold);
RcppExport SEXP _dap_dap_main(SEXP XSEXP, SEXP ySEXP, SEXP matrixSEXP, SEXP prior_weightsSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP, SEXP phi2_matSEXP, SEXP exclusiveSEXP, SEXP pir_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_weights(prior_weightsSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi2_mat(phi2_matSEXP);
    Rcpp::traits::input_parameter< bool >::type exclusive(exclusiveSEXP);
    Rcpp::traits::input_parameter< double >::type pir_threshold(pir_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(dap_main(X, y, matrix, prior_weights, r2_threshold, coverage, phi2_mat, exclusive, pir_threshold));
    return rcpp_result_gen;
END_RCPP
}
// pir
List pir(NumericMatrix mat, double threshold);
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
List get_sc(const NumericMatrix& X, const NumericMatrix& combo, const NumericMatrix& single, const NumericVector& posterior_prob, const CharacterVector& col_names, double threshold, double r2_threshold, double coverage);
RcppExport SEXP _dap_get_sc(SEXP XSEXP, SEXP comboSEXP, SEXP singleSEXP, SEXP posterior_probSEXP, SEXP col_namesSEXP, SEXP thresholdSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type combo(comboSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type single(singleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type posterior_prob(posterior_probSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type col_names(col_namesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sc(X, combo, single, posterior_prob, col_names, threshold, r2_threshold, coverage));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_compute_log10_prior", (DL_FUNC) &_dap_compute_log10_prior, 2},
    {"_dap_compute_log10_posterior", (DL_FUNC) &_dap_compute_log10_posterior, 6},
    {"_dap_dap_main", (DL_FUNC) &_dap_dap_main, 9},
    {"_dap_pir", (DL_FUNC) &_dap_pir, 2},
    {"_dap_get_sc", (DL_FUNC) &_dap_get_sc, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
