// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_log10_posterior
List compute_log10_posterior(const NumericMatrix& X, const NumericVector& y, const std::vector<std::vector<int>>& cmfg_matrix, const NumericVector& pi_vec, const NumericMatrix& phi2_mat);
RcppExport SEXP _dap_compute_log10_posterior(SEXP XSEXP, SEXP ySEXP, SEXP cmfg_matrixSEXP, SEXP pi_vecSEXP, SEXP phi2_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type cmfg_matrix(cmfg_matrixSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type phi2_mat(phi2_matSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_posterior(X, y, cmfg_matrix, pi_vec, phi2_mat));
    return rcpp_result_gen;
END_RCPP
}
// dap_main
List dap_main(NumericMatrix X, NumericVector y, NumericMatrix matrix, double pir_threshold, NumericVector prior_weights, NumericMatrix phi2_mat, double r2_threshold, double coverage, bool exclusive, CharacterVector snp_names);
RcppExport SEXP _dap_dap_main(SEXP XSEXP, SEXP ySEXP, SEXP matrixSEXP, SEXP pir_thresholdSEXP, SEXP prior_weightsSEXP, SEXP phi2_matSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP, SEXP exclusiveSEXP, SEXP snp_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< double >::type pir_threshold(pir_thresholdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_weights(prior_weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi2_mat(phi2_matSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< bool >::type exclusive(exclusiveSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type snp_names(snp_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(dap_main(X, y, matrix, pir_threshold, prior_weights, phi2_mat, r2_threshold, coverage, exclusive, snp_names));
    return rcpp_result_gen;
END_RCPP
}
// get_sc
List get_sc(const NumericMatrix& X, const NumericMatrix& effect_pip, const CharacterVector& snp_names, double r2_threshold, double coverage);
RcppExport SEXP _dap_get_sc(SEXP XSEXP, SEXP effect_pipSEXP, SEXP snp_namesSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type effect_pip(effect_pipSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type snp_names(snp_namesSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sc(X, effect_pip, snp_names, r2_threshold, coverage));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_compute_log10_posterior", (DL_FUNC) &_dap_compute_log10_posterior, 5},
    {"_dap_dap_main", (DL_FUNC) &_dap_dap_main, 10},
    {"_dap_get_sc", (DL_FUNC) &_dap_get_sc, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
