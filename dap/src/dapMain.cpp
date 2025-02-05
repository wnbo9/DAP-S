#include <Rcpp.h>
#include <iostream>     // for cerr, endl
#include <fstream>      // for ofstream
#include <vector>       // for vector
#include <string>       // for string
#include <utility>      // for
#include <vector>

using namespace Rcpp;
using namespace std;

vector<vector<int>> pir(const vector<vector<double>>& mat, double pir_threshold);
List compute_log10_posterior(const NumericMatrix& X, const NumericVector& y, const std::vector<std::vector<int>>& cmfg_matrix, const NumericVector& pi_vec, const NumericMatrix& phi2_mat);
List get_sc(const NumericMatrix& X, const NumericMatrix& effect_pip, const CharacterVector& snp_names, double r2_threshold, double coverage);

//' Implementation of DAP-S algorithm in C++ with default SuSiE settings
//' 
//' Performs fine mapping using DAP-S algorithm
//'
//' @param X Genotype matrix
//' @param y Phenotype vector
//' @param matrix Proposal density matrix from SuSiE
//' @param pir_threshold Threshold for PIR
//' @param prior_weights Vector of prior probabilities
//' @param phi2_mat Matrix of scaled prior effect size variances
//' @param r2_threshold Threshold for LD
//' @param coverage Coverage for credible set
//' @param exclusive If TRUE. enforce mutually exclusive clusters
//' @param snp_names SNP names
//'
//' @return A list containing:
//' \itemize{
//'   \item model_config - Model configurations
//'   \item posterior_prob - Posterior probabilities
//'   \item log10_posterior_score - Log10 posterior scores
//'   \item log10_nc - Log10 normalizing constant
//'   \item pip - Posterior inclusion probabilities
//'   \item signal_cluster - Signal clusters
//' }
//' @export
// [[Rcpp::export]]
List dap_main(NumericMatrix X, 
              NumericVector y,
              NumericMatrix matrix,
              double pir_threshold,
              NumericVector prior_weights,
              NumericMatrix phi2_mat,
              double r2_threshold,
              double coverage,
              bool exclusive,
              CharacterVector snp_names) {
    
    int p = X.ncol(); // 5000
    int L = matrix.ncol(); // 3

    vector<vector<double>> mat(p+1, vector<double>(L));
    for (int i = 0; i < p+1; i++) {
        for (int j = 0; j < L; j++) {
            mat[i][j] = matrix(i, j);
        }
    }

    // PIR
    Rcout << "---Pseudo Importance Resampling...\n";
    Rcout << "---PIR threshold: " << pir_threshold << endl;

    vector<vector<int>> combo_matrix = pir(mat, pir_threshold);
    int m_size = combo_matrix.size();
    NumericMatrix combo(m_size, L);
    for(int i = 0; i < m_size; i++) {
        for(size_t j = 0; j < L; j++) {
            combo(i, j) = combo_matrix[i][j];
        }
    }

    // Compute log10 posterior
    Rcpp::Rcout << "---Calculating posterior of " << m_size << " model configurations...\n";
    List scores = compute_log10_posterior(X, y, combo_matrix, prior_weights, phi2_mat);
    NumericVector log10_BF = scores["log10_BF"];
    NumericVector log10_prior = scores["log10_prior"];
    NumericVector log10_posterior_score = scores["log10_posterior_score"];

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

    // Print model configurations
    CharacterVector model_config(m_size);
    for (int i = 0; i < m_size; i++) {
        vector<int>& row = combo_matrix[i];
        string config = "";
        bool all_p = true;
    
        for (int j = 0; j < L; j++) {
            if (row[j] != p) {
                all_p = false;
                if (config != "") {
                    config += "+";
                }
                config += string(snp_names[row[j]]);
            }
        }
    
        if (all_p) {
            model_config[i] = "NULL";
        } else {
            model_config[i] = config;
        }
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


    // Construct signal clusters
    Rcout << "---Constructing " << (exclusive ? "mutually exclusive " : "") 
          << (coverage > 1 ? "signal clusters" : (coverage > 0 ? to_string(int(coverage * 100)) + "% credible sets" : "? credible sets")) 
          << endl;
    List sc_results = get_sc(X, effect_pip, snp_names, r2_threshold, coverage);

    // Return simplified results
    return List::create(
        Named("model_config") = model_config,
        Named("posterior_prob") = posterior_prob,
        Named("log10_posterior_score") = log10_posterior_score,
        Named("log10_BF") = log10_BF,
        Named("log10_prior") = log10_prior,
        Named("log10_nc") = log_nc,
        Named("pip") = marginal_pip,
        Named("effect_pip") = effect_pip,
        Named("signal_cluster") = sc_results
    );
}
