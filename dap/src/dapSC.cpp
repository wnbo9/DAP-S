#include <Rcpp.h>
#include <gsl/gsl_statistics.h>
#include <vector>
#include <set>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(gsl)]]

//' Compute R-squared between two vectors
//' 
//' @param X NumericMatrix containing the vectors
//' @param i Index of first vector
//' @param j Index of second vector
//' @return Double containing R-squared value
//' @noRd
double compute_r2(const NumericMatrix& X, int i, int j) {
    double r2;
    int n = X.nrow();
    const double* gi = &X(0, i);
    const double* gj = &X(0, j);
    
    r2 = pow(gsl_stats_correlation(gi, 1, gj, 1, n), 2);
    return r2;
}


//' Get signal clusters or credible sets at given coverage level
//' 
//' This function identifies clusters of correlated signals based on R-squared values
//' and model posterior probabilities.
//' 
//' @param X NumericMatrix containing the raw data
//' @param mat NumericMatrix containing the alpha matrix
//' @param cmfg_mat NumericMatrix containing the CMFG matrix
//' @param posterior_probs NumericVector of posterior probabilities
//' @param col_names CharacterVector of column names
//' @param r2_threshold Double specifying the R-squared threshold for correlation
//' @param coverage Double specifying the coverage threshold
//' @return List containing:
//'   \item{clusters}{List of character vectors containing cluster memberships}
//'   \item{spip}{Numeric vector of signal posterior inclusion probabilities}
//'   \item{sizes}{Integer vector of cluster sizes}
//'   \item{r2_threshold}{Double containing the R-squared threshold used}
//'   \item{coverage}{Double containing the coverage threshold used}
//' @export
// [[Rcpp::export]]
List get_sc(const NumericMatrix& X,
            const NumericMatrix& mat,
            const NumericMatrix& cmfg_mat,
            const NumericVector& posterior_probs,
            const CharacterVector& col_names,
            double r2_threshold,
            double coverage) {
    
    int p = mat.nrow();
    int L = mat.ncol();
    int m = cmfg_mat.nrow();
    std::vector<std::vector<int>> clusters;
    std::vector<double> mps_values;


    // Sort each column in descending order
    vector<vector<pair<double, int>>> sortedMat(L);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < p; j++) {
            sortedMat[i].emplace_back(mat(j, i), j);
        }
        sort(sortedMat[i].begin(), sortedMat[i].end(), greater<pair<double, int>>());
    }


    // For each column in alpha matrix
    for (int col = 0; col < L; col++) {
        // Skip if the null SNP is top in the column
        int top_snp = sortedMat[col][0].second;
        if(top_snp == p - 1) continue;

        // Start with the top SNP
        std::vector<int> current_cluster = {top_snp};

        // Initialize vector to track models with no SNPs in it
        std::vector<int> models_without_snps;
        double sum_prob_without = 0.0;

        // Initially find models without top_snp
        for(int i = 0; i < m; i++) {
            if(cmfg_mat(i, top_snp) == 0) {
                models_without_snps.push_back(i);
                sum_prob_without += posterior_probs[i];
            }
        }
        double mps = 1.0 - sum_prob_without;

        // Add SNP
        for (size_t j = 1; j < p; j++){
            int next_snp = sortedMat[col][j].second;
            if (next_snp == p - 1) continue;

            bool add_to_cluster = true;

            // Check R2 with all SNPs in current cluster
            for (int cluster_snp : current_cluster) {
                if (compute_r2(X, cluster_snp, next_snp) < r2_threshold) {
                    add_to_cluster = false;
                    break;
                }
            }

            // Calculate MPS for current cluster
            if (add_to_cluster){
                current_cluster.push_back(next_snp);

                // Calculate MPS for current cluster
                for (size_t idx = 0; idx < models_without_snps.size();) {
                    int model_idx = models_without_snps[idx];
                    if (cmfg_mat(model_idx, next_snp) == 1) {
                        sum_prob_without -= posterior_probs[model_idx];
                        models_without_snps[idx] = models_without_snps.back();
                        models_without_snps.pop_back();
                    } else {
                        idx++;
                    }
                }
                mps = 1 - sum_prob_without;

                // If MPS exceeds coverage threshold, add cluster and break
                if (mps >= coverage) {
                    break;
                }
            }
        }

        // After going through all SNPs, add the final cluster
        if(!current_cluster.empty()) {
            clusters.push_back(current_cluster);
            mps_values.push_back(mps);
        }
    }

    // Only include non-empty clusters in results
    std::vector<std::vector<int>> non_empty_clusters;
    std::vector<double> non_empty_mps;
    for(size_t i = 0; i < clusters.size(); i++) {
        if(!clusters[i].empty()) {
            non_empty_clusters.push_back(clusters[i]);
            non_empty_mps.push_back(mps_values[i]);
        }
    }

    // Convert results to R objects
    int n_clusters = non_empty_clusters.size();
    List r_clusters(n_clusters);
    NumericVector r_mps(n_clusters);
    IntegerVector cluster_sizes(n_clusters);
    
    for(int i = 0; i < n_clusters; i++) {
        CharacterVector cluster_names(non_empty_clusters[i].size());
        for(size_t j = 0; j < non_empty_clusters[i].size(); j++) {
            cluster_names[j] = col_names[non_empty_clusters[i][j]];
        }
        
        r_clusters[i] = cluster_names;
        r_mps[i] = non_empty_mps[i];
        cluster_sizes[i] = non_empty_clusters[i].size();
    }
    
    return List::create(
        Named("clusters") = r_clusters,
        Named("spip") = r_mps,
        Named("sizes") = cluster_sizes,
        Named("r2_threshold") = r2_threshold,
        Named("coverage") = coverage
    );
}




