#include <Rcpp.h>
#include <iostream>     // for cerr, endl
#include <fstream>      // for ofstream
#include <vector>       // for vector
#include <string>       // for string
#include <utility>      // for
#include <vector>

using namespace Rcpp;
using namespace std;

double compute_log10_prior(const std::vector<int> &mcfg, NumericVector pi_vec);
List get_sc(const NumericMatrix& X, const NumericMatrix& effect_pip, const CharacterVector& snp_names, double r2_threshold, double coverage);

//' Update DAP-S results
//' @param X Genotype matrix
//' @param dap_result DAP-S results
//' @param prior_weights Vector of prior probabilities
//' @param r2_threshold Threshold for LD
//' @param coverage Coverage for credible set
//'
//' @return A list containing:
//' \itemize{
//'   \item posterior_prob - Posterior probabilities
//'   \item log10_posterior_score - Log10 posterior scores
//'   \item log10_nc - Log10 normalizing constant
//'   \item pip - Posterior inclusion probabilities
//'   \item signal_cluster - Signal clusters
//' }
// [[Rcpp::export]]
List dap_update_main(NumericMatrix X, 
                     List dap_result, 
                     NumericVector prior_weights,
                     double r2_threshold,
                     double coverage) {

    List models = dap_result["models"];
    NumericVector log10_BF = models["log10_BF"];
    vector<vector<int>> combo_matrix = dap_result["model_combo"];
    int m_size = combo_matrix.size();
    int L = combo_matrix[0].size();
    CharacterVector snp_names = dap_result["snp_names"];
    List params = dap_result["params"];
    int p = snp_names.size();
    bool overlapping = params["overlapping"];

    NumericVector log10_prior(m_size);
    NumericVector log10_posterior_score(m_size);

    for (int i = 0; i < m_size; i++) {
        log10_prior[i] = compute_log10_prior(combo_matrix[i], prior_weights);
        log10_posterior_score[i] = log10_BF[i] + log10_prior[i];
    }

    double max_log_posterior = *max_element(log10_posterior_score.begin(), log10_posterior_score.end());
    double sum_exp = 0.0;
    NumericVector posterior_prob(m_size);
    for(int i = 0; i < m_size; i++) {
        sum_exp += pow(10.0, log10_posterior_score[i] - max_log_posterior);
    }
    double log_nc = max_log_posterior + log10(sum_exp);

    // Model posterior probabilities
    for(int i = 0; i < m_size; i++) {
        posterior_prob[i] = pow(10.0, log10_posterior_score[i] - log_nc);
    }

    // Calculate effect PIP and marginal PIP
    NumericMatrix effect_pip(p, L);
    NumericVector marginal_pip(p);
    for(int i = 0; i < m_size; i++) { // For each combination
        for(int j = 0; j < L; j++) {
            if(combo_matrix[i][j] < p) {
                effect_pip(combo_matrix[i][j], j) += posterior_prob[i];
                marginal_pip[combo_matrix[i][j]] += posterior_prob[i];
            }
        }
    }

    // Clamp values between 0 and 1
    for(int i = 0; i < p; i++) {
        for (int j = 0; j < L; j++) {
            effect_pip(i, j) = max(0.0, min(1.0, effect_pip(i, j)));
        }
        marginal_pip[i] = max(0.0, min(1.0, marginal_pip[i]));
    }
    rownames(effect_pip) = snp_names;
    marginal_pip.names() = snp_names;

    // Calculate TWAS weights
    bool twas_weight = params["twas_weight"];
    NumericVector twas_weights(p);
    NumericMatrix reg_weights = dap_result["reg_weights"];
    if (twas_weight) {
        const double* post_ptr = posterior_prob.begin();
        for(int j = 0; j < p; j++) {
            double sum = 0.0;
            for(int i = 0; i < m_size; i++) {
                sum += post_ptr[i] * reg_weights(i,j);
            }
            twas_weights[j] = sum;
        }
    }

    // Construct signal clusters
    Rcout << "---Constructing " << (overlapping ? "overlapping " : "") 
          << (coverage > 1 ? "signal clusters" : (coverage > 0 ? to_string(int(coverage * 100)) + "% credible sets" : "? credible sets")) 
          << endl;
    List sc_results = get_sc(X, effect_pip, snp_names, r2_threshold, coverage);

    // Return simplified results
    return List::create(
        Named("posterior_prob") = posterior_prob, //
        Named("log10_posterior_score") = log10_posterior_score, //
        Named("log10_BF") = log10_BF,
        Named("log10_prior") = log10_prior, //
        Named("log10_nc") = log_nc, //
        Named("pip") = marginal_pip, //
        Named("effect_pip") = effect_pip, //
        Named("signal_cluster") = sc_results, //
        Named("twas_weights") = twas_weights //
    );
}