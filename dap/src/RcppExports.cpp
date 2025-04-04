// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_log10_posterior
List compute_log10_posterior(const std::vector<std::vector<int>>& cmfg_matrix, const NumericVector& pi_vec, const NumericMatrix& phi2_mat, int ss, SEXP X_input, SEXP y_input, SEXP XtX_input, SEXP Xty_input, SEXP yty_input, SEXP n_input, SEXP V_input, SEXP Dsq_input, SEXP var_input, SEXP XtOmegay_input, bool twas_weight);
RcppExport SEXP _dap_compute_log10_posterior(SEXP cmfg_matrixSEXP, SEXP pi_vecSEXP, SEXP phi2_matSEXP, SEXP ssSEXP, SEXP X_inputSEXP, SEXP y_inputSEXP, SEXP XtX_inputSEXP, SEXP Xty_inputSEXP, SEXP yty_inputSEXP, SEXP n_inputSEXP, SEXP V_inputSEXP, SEXP Dsq_inputSEXP, SEXP var_inputSEXP, SEXP XtOmegay_inputSEXP, SEXP twas_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type cmfg_matrix(cmfg_matrixSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type phi2_mat(phi2_matSEXP);
    Rcpp::traits::input_parameter< int >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type X_input(X_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_input(y_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type XtX_input(XtX_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Xty_input(Xty_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yty_input(yty_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_input(n_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type V_input(V_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Dsq_input(Dsq_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type var_input(var_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type XtOmegay_input(XtOmegay_inputSEXP);
    Rcpp::traits::input_parameter< bool >::type twas_weight(twas_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10_posterior(cmfg_matrix, pi_vec, phi2_mat, ss, X_input, y_input, XtX_input, Xty_input, yty_input, n_input, V_input, Dsq_input, var_input, XtOmegay_input, twas_weight));
    return rcpp_result_gen;
END_RCPP
}
// dap_main
List dap_main(NumericMatrix matrix, double pir_threshold, NumericVector prior_weights, NumericMatrix phi2_mat, double r2_threshold, double coverage, bool overlapping, bool twas_weight, CharacterVector snp_names, int ss, SEXP X_input, SEXP y_input, SEXP XtX_input, SEXP Xty_input, SEXP yty_input, SEXP n_input, SEXP V_input, SEXP Dsq_input, SEXP var_input, SEXP XtOmegay_input);
RcppExport SEXP _dap_dap_main(SEXP matrixSEXP, SEXP pir_thresholdSEXP, SEXP prior_weightsSEXP, SEXP phi2_matSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP, SEXP overlappingSEXP, SEXP twas_weightSEXP, SEXP snp_namesSEXP, SEXP ssSEXP, SEXP X_inputSEXP, SEXP y_inputSEXP, SEXP XtX_inputSEXP, SEXP Xty_inputSEXP, SEXP yty_inputSEXP, SEXP n_inputSEXP, SEXP V_inputSEXP, SEXP Dsq_inputSEXP, SEXP var_inputSEXP, SEXP XtOmegay_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< double >::type pir_threshold(pir_thresholdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_weights(prior_weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi2_mat(phi2_matSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< bool >::type overlapping(overlappingSEXP);
    Rcpp::traits::input_parameter< bool >::type twas_weight(twas_weightSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type snp_names(snp_namesSEXP);
    Rcpp::traits::input_parameter< int >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type X_input(X_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y_input(y_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type XtX_input(XtX_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Xty_input(Xty_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yty_input(yty_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_input(n_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type V_input(V_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Dsq_input(Dsq_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type var_input(var_inputSEXP);
    Rcpp::traits::input_parameter< SEXP >::type XtOmegay_input(XtOmegay_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(dap_main(matrix, pir_threshold, prior_weights, phi2_mat, r2_threshold, coverage, overlapping, twas_weight, snp_names, ss, X_input, y_input, XtX_input, Xty_input, yty_input, n_input, V_input, Dsq_input, var_input, XtOmegay_input));
    return rcpp_result_gen;
END_RCPP
}
// dap_update_main
List dap_update_main(NumericMatrix X, List dap_result, NumericVector prior_weights, double r2_threshold, double coverage);
RcppExport SEXP _dap_dap_update_main(SEXP XSEXP, SEXP dap_resultSEXP, SEXP prior_weightsSEXP, SEXP r2_thresholdSEXP, SEXP coverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type dap_result(dap_resultSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_weights(prior_weightsSEXP);
    Rcpp::traits::input_parameter< double >::type r2_threshold(r2_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type coverage(coverageSEXP);
    rcpp_result_gen = Rcpp::wrap(dap_update_main(X, dap_result, prior_weights, r2_threshold, coverage));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_compute_log10_posterior", (DL_FUNC) &_dap_compute_log10_posterior, 15},
    {"_dap_dap_main", (DL_FUNC) &_dap_dap_main, 20},
    {"_dap_dap_update_main", (DL_FUNC) &_dap_dap_update_main, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
