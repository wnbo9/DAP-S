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
List compute_log10_posterior(const NumericMatrix& X, const NumericVector& y, const std::vector<std::vector<int>>& cmfg_matrix, const NumericVector& pi_vec, const NumericMatrix& phi2_mat, bool twas_weight);
List get_sc(const NumericMatrix& X, const NumericMatrix& effect_pip, const CharacterVector& snp_names, double r2_threshold, double coverage);

//' Implementation of DAP-S algorithm in C++ with default SuSiE settings
//' @param X Genotype matrix
//' @param y Phenotype vector
//' @param matrix Proposal density matrix from SuSiE
//' @param pir_threshold Threshold for PIR
//' @param prior_weights Vector of prior probabilities
//' @param phi2_mat Matrix of scaled prior effect size variances
//' @param r2_threshold Threshold for LD
//' @param coverage Coverage for credible set
//' @param overlapping If TRUE. enforce overlapping clusters
//' @param twas_weight If TRUE, return TWAS weights
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
// [[Rcpp::export]]
List dap_main(NumericMatrix X, 
              NumericVector y,
              NumericMatrix matrix,
              double pir_threshold,
              NumericVector prior_weights,
              NumericMatrix phi2_mat,
              double r2_threshold,
              double coverage,
              bool overlapping,
              bool twas_weight,
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

    // Compute log10 posterior
    Rcpp::Rcout << "---Calculating posterior of " << m_size << " model configurations...\n";
    List scores = compute_log10_posterior(X, y, combo_matrix, prior_weights, phi2_mat, twas_weight);
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
        bool first = true; // first is null

        for (int snp : row) {
            if (snp != p) {
                if (!first) {
                    config += "+";
                }
                config += string(snp_names[snp]);
                first = false;
            }
        }
        model_config[i] = first ? "NULL" : config;
    }

    // Calculate effect PIP and marginal PIP
    NumericMatrix effect_pip(p, L);
    NumericVector marginal_pip(p);

    // First create a vector storing the SNP at each effect position
    vector<set<int>> effects_snp(L); // initialize with p (null)
    for (int i = 0; i < m_size; i++) {
        vector<int>& model = combo_matrix[i];
        if (model.size() > 1) {  // only look at multi-SNP models
            for (int pos = 0; pos < model.size(); pos++) {
                int snp = model[pos];
                if (snp != p) {
                    effects_snp[pos].insert(snp);
                }
            }
        }
    }

    for(int i = 0; i < m_size; i++) { // For each combination
        vector<int>& model = combo_matrix[i];

        if (model.size() == 1) {
            // For 1-snp models, add to all columns of effect PIP and update marginal PIP
            if (model[0] != p) {
                marginal_pip[model[0]] += posterior_prob[i];
                for (int l = 0; l < L; l++) {
                    if (effects_snp[l].find(model[0]) != effects_snp[l].end()) {
                        effect_pip(model[0], l) += posterior_prob[i];
                    }
                }
            }
        } else {
            // For multi-snp models
            for (int pos = 0; pos < model.size(); pos++) {
                int snp = model[pos];
                if (snp != p) {
                    marginal_pip[snp] += posterior_prob[i];
                    effect_pip(snp, pos) += posterior_prob[i];
                }
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
    if (L == 1) {
        for (int i = 0; i < p; i++) {
            effect_pip(i, 0) = marginal_pip[i];
        }
    }
    rownames(effect_pip) = snp_names;

    // Calculate TWAS weights
    NumericVector twas_weights(p);
    NumericMatrix reg_weights = scores["reg_weights"];
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
        Named("model_config") = model_config,
        Named("model_combos") = combo_matrix,
        Named("posterior_prob") = posterior_prob,
        Named("log10_posterior_score") = log10_posterior_score,
        Named("log10_BF") = log10_BF,
        Named("log10_prior") = log10_prior,
        Named("log10_nc") = log_nc,
        Named("pip") = marginal_pip,
        Named("effect_pip") = effect_pip,
        Named("signal_cluster") = sc_results,
        Named("reg_weights") = reg_weights,
        Named("twas_weights") = twas_weights
    );
}