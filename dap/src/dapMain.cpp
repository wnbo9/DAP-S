#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

List pir(NumericMatrix mat, double threshold);
List compute_log10_posterior(NumericMatrix X, NumericVector y, NumericMatrix binaryCombinations, NumericVector pi_vec, NumericVector phi2_vec);
List get_sc(const NumericMatrix& X, const NumericMatrix& mat, const NumericMatrix& cmfg_mat, const NumericVector& posterior_prob, const CharacterVector& col_names, double threshold, double r2_threshold, double coverage);

//' Implementation of DAP-PIR algorithm in C++
//' 
//' Performs fine mapping using DAP-PIR algorithm
//'
//' @param X Genotype matrix
//' @param y Phenotype vector
//' @param matrix Proposal density matrix from SuSiE
//' @param threshold Threshold for proposal density
//' @param prior_weights Vector of prior probabilities
//' @param phi2_vec Vector of scaled prior effect size variances
//' @param r2_threshold Threshold for LD
//' @param coverage Coverage for credible set
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
              NumericMatrix matrix,
              double threshold,
              NumericVector prior_weights,
              NumericVector phi2_vec,
              double r2_threshold,
              double coverage) {

    // Get the number of rows and columns
    int p = X.ncol();
    CharacterVector col_names = colnames(X);

    // PIR
    Rcpp::Rcout << "Pseudo Importance Resamping...\n";
    List results = pir(matrix, threshold);
    List positionList = results["position_elements"];
    NumericMatrix cmfg_mat = results["combinations"];

    int m_size = cmfg_mat.nrow();

    // Compute log10 posterior
    Rcpp::Rcout << "Calculating posterior of " << m_size << " model configurations...\n";
    List scores = compute_log10_posterior(X, y, cmfg_mat, prior_weights, phi2_vec);
    NumericVector log10_BF = scores["log10_BF"];
    NumericVector log10_prior = scores["log10_prior"];
    NumericVector log10_posterior_score = scores["log10_posterior_score"];

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

    // Calculate PIP
    NumericVector pip(p);
    for(int i = 0; i < m_size; i++) {
        for(int j = 0; j < p; j++) {
            if(cmfg_mat(i, j) == 1) {
                pip[j] += posterior_prob[i];
            }
        }
    }
    // Ensure PIPs are bounded between 0 and 1 (numerical stability)
    for(int j = 0; j < p; j++) {
        if(pip[j] > 1.0) {
            pip[j] = 1.0;
        } else if(pip[j] < 0.0) {
            pip[j] = 0.0;
        }
    }
    pip.names() = col_names;


    // Print model configurations
    CharacterVector model_config(m_size);
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
            model_config[i] = "NULL";
        } else {
            model_config[i] = config;
        }
    }
    

    // Construct signal clusters
    List sc_results = get_sc(X, matrix, cmfg_mat, posterior_prob, col_names, threshold, r2_threshold, coverage);

    // Return simplified results
    return List::create(
        Named("model_config") = model_config,
        Named("posterior_prob") = posterior_prob,
        Named("log10_posterior_score") = log10_posterior_score,
        Named("log10_BF") = log10_BF,
        Named("log10_prior") = log10_prior,
        Named("log10_nc") = log_nc,
        Named("pip") = pip,
        Named("signal_cluster") = sc_results,
        Named("element_cluster") = positionList
    );
}
