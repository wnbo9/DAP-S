#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

List pir(NumericMatrix mat, double pir_threshold);
List compute_log10_posterior(NumericMatrix X, NumericVector y, NumericMatrix matrixCombinations, NumericMatrix matrixSingle, NumericVector pi_vec, NumericMatrix phi2_mat);
List get_sc(const NumericMatrix& X, const NumericMatrix& combo, const NumericMatrix& single, const NumericVector& posterior_prob, const CharacterVector& col_names, double threshold, double r2_threshold, double coverage);

//' Implementation of DAP-S algorithm in C++ with default SuSiE settings
//' 
//' Performs fine mapping using DAP-S algorithm
//'
//' @param X Genotype matrix
//' @param y Phenotype vector
//' @param matrix Proposal density matrix from SuSiE
//' @param prior_weights Vector of prior probabilities
//' @param r2_threshold Threshold for LD
//' @param coverage Coverage for credible set
//' @param phi2_mat Matrix of scaled prior effect size variances
//' @param exclusive If TRUE. enforce mutually exclusive clusters
//' @param pir_threshold Threshold for PIR
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
              NumericVector prior_weights,
              double r2_threshold,
              double coverage,
              NumericMatrix phi2_mat,
              bool exclusive,
              double pir_threshold) {
    
    Rcout << "---PIR threshold: " << pir_threshold << std::endl;

    // Get the number of rows and columns
    int p = X.ncol(); // 5000
    CharacterVector col_names = colnames(X);

    // PIR
    Rcpp::Rcout << "---Pseudo Importance Resamping...\n";
    List results = pir(matrix, pir_threshold);
    NumericMatrix combo = results["combo"];
    NumericMatrix single = results["single"];

    int L = combo.ncol();
    int m1_size = combo.nrow();
    int m2_size = single.nrow();
    int m_size = m1_size + m2_size;
    Rcout << "---Number of combo configurations: " << m1_size << std::endl;
    Rcout << "---Number of single configurations: " << m2_size << std::endl;


    // Compute log10 posterior
    Rcpp::Rcout << "---Calculating posterior of " << m_size << " model configurations...\n";
    List scores = compute_log10_posterior(X, y, combo, single, prior_weights, phi2_mat);
    NumericVector log10_BF = scores["log10_BF"];
    NumericVector log10_prior = scores["log10_prior"];
    NumericVector log10_posterior_score = scores["log10_posterior_score"];
    

    if (log10_BF.length() != m_size || log10_prior.length() != m_size || log10_posterior_score.length() != m_size) {
    stop("Mismatch in vector sizes. Expected size: %d, Got: %d, %d, %d", 
         m_size, log10_BF.length(), log10_prior.length(), log10_posterior_score.length());
    }


    // Calculate normalizing constant
    double max_log_posterior = *std::max_element(log10_posterior_score.begin(), log10_posterior_score.end());
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
        NumericVector row;
        if (i < m1_size) {
            row = combo(i, _);
        } else {
            row = single(i - m1_size, _);
        }
        std::string config = "";
        bool all_p = true;
    
        for (int j = 0; j < row.size(); j++) {
            if (row[j] != p) {
                all_p = false;
                if (config != "") {
                    config += "+";
                }
                config += std::string(col_names[row[j]]);
            }
        }
    
        if (all_p) {
            model_config[i] = "NULL";
        } else {
            model_config[i] = config;
        }
    }

    // Calculate marginal PIP
    NumericVector marginal_pip(p);
    for(int i = 0; i < m1_size; i++) { // For each combination
        for(int j = 0; j < L; j++) {
            if(combo(i, j) < p) marginal_pip[combo(i, j)] += posterior_prob[i];
        }
    }
    for(int i = 0; i < m2_size; i++) { // For each single SNP
        if (single(i, 0) < p) marginal_pip[single(i, 0)] += posterior_prob[m1_size + i];
    }
    // Ensure PIPs are bounded between 0 and 1 (numerical stability)
    for(int j = 0; j < p; j++) {
        if(marginal_pip[j] > 1.0) {
            marginal_pip[j] = 1.0;
        } else if(marginal_pip[j] < 0.0) {
            marginal_pip[j] = 0.0;
        }
    }
    marginal_pip.names() = col_names;


    // Construct signal clusters
    Rcout << "---Constructing " << (exclusive ? "mutually exclusive " : "") 
          << (coverage > 1 ? "signal clusters" : (coverage > 0 ? std::to_string(int(coverage * 100)) + "% credible sets" : "? credible sets")) 
          << std::endl;
    List sc_results = get_sc(X, combo, single, posterior_prob, col_names, pir_threshold, r2_threshold, coverage);

    // Return simplified results
    return List::create(
        Named("model_config") = model_config,
        Named("posterior_prob") = posterior_prob,
        Named("log10_posterior_score") = log10_posterior_score,
        Named("log10_BF") = log10_BF,
        Named("log10_prior") = log10_prior,
        Named("log10_nc") = log_nc,
        Named("pip") = marginal_pip,
        Named("signal_cluster") = sc_results,
        Named("combo") = combo,
        Named("single") = single
    );
}
