#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

NumericMatrix pir(NumericMatrix mat, double threshold);
NumericVector compute_log10_posterior(NumericMatrix X, NumericVector y,
                                    NumericMatrix binaryCombinations, 
                                    NumericVector pi_vec, NumericVector phi2_vec);

//' Main function for DAP-PIR algorithm
//' 
//' Performs fine mapping using DAP-PIR algorithm
//'
//' @param X Genotype matrix
//' @param y Phenotype vector
//' @param L Number of causal variants
//' @param matrix Proposal density matrix from SuSiE
//' @param threshold Threshold for proposal density
//' @param prior_weights Vector of prior probabilities
//' @param phi2_vec Vector of scaled prior effect size variances
//'
//' @return A list containing:
//' \itemize{
//'   \item model_config - Model configurations
//'   \item posterior_prob - Posterior probabilities
//'   \item log10_posterior_score - Log10 posterior scores
//'   \item log10_nc - Log10 normalizing constant
//'   \item pip - Posterior inclusion probabilities
//' }
//' @export
// [[Rcpp::export]]
List dap_main(NumericMatrix X, 
              NumericVector y,
              int L,
              NumericMatrix matrix,
              double threshold,
              NumericVector prior_weights,
              NumericVector phi2_vec) {

    // Get the number of rows and columns
    int n = X.nrow();
    int p = X.ncol();
    CharacterVector col_names = colnames(X);

    // PIR
    Rcpp::Rcout << "Pseudo Importance Resamping...\n";
    NumericMatrix cmfg_mat = pir(matrix, threshold);
    int m_size = cmfg_mat.nrow();

    // Compute log10 posterior
    Rcpp::Rcout << "Calculating posterior of " << m_size << " model configurationss...\n";
    NumericVector log10_posterior = compute_log10_posterior(X, y, cmfg_mat, prior_weights, phi2_vec);

    // Calculate normalizing constant
    double max_log_posterior = *std::max_element(log10_posterior.begin(), log10_posterior.end());
    double sum_exp = 0.0;
    NumericVector posterior_probs(m_size);
    for(int i = 0; i < m_size; i++) {
        sum_exp += pow(10.0, log10_posterior[i] - max_log_posterior);
    }
    double log_nc = max_log_posterior + log10(sum_exp);

    // Model posterior probabilities
    for(int i = 0; i < m_size; i++) {
        posterior_probs[i] = pow(10.0, log10_posterior[i] - log_nc);
    }

    // Calculate PIP
    NumericVector pip(p);
    for(int i = 0; i < m_size; i++) {
        for(int j = 0; j < p; j++) {
            if(cmfg_mat(i, j) == 1) {
                pip[j] += posterior_probs[i];
            }
        }
    }
    pip.names() = col_names;

    // Print model configurations
    CharacterVector model_configs(m_size);
    for(int i = 0; i < m_size; i++) {
        string config = "";
        bool first = true;
        bool all_zero = true;
        
        for(int j = 0; j < p; j++) {
            if(cmfg_mat(i, j) == 1) {
                all_zero = false;
                if(!first) {
                    config += " + ";
                }
                config += as<string>(col_names[j]);
                first = false;
            }
        }
        
        if(all_zero) {
            model_configs[i] = "NULL";
        } else {
            model_configs[i] = config;
        }
    }
    

    // Construct signal clusters
    // Compute signal level PIP

    // Construct credible sets




    // Return simplified results
    return List::create(
        Named("model_config") = model_configs,
        Named("posterior_prob") = posterior_probs,
        Named("log10_posterior_score") = log10_posterior,
        Named("log10_nc") = log_nc,
        Named("pip") = pip
    );
}

