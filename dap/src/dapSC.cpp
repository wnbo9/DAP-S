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

//' Compute median of a sorted vector
//' @name median
//' @param v Vector of sorted values
//' @return Median value
//' @noRd
double median(vector<double>& v) {
    size_t n = v.size();
    if(n % 2 == 0) {
        return (v[n/2 - 1] + v[n/2]) / 2;
    }
    return v[n/2];
}

//' Get signal clusters or credible sets at given coverage level
//' 
//' This function identifies clusters of correlated signals based on R-squared values
//' and model posterior probabilities.
//' 
//' @param X NumericMatrix containing the raw data
//' @param mat NumericMatrix containing the alpha matrix
//' @param cmfg_mat NumericMatrix containing the CMFG matrix
//' @param posterior_prob NumericVector of posterior probabilities
//' @param col_names CharacterVector of column names
//' @param threshold Double specifying the threshold for model proposal density
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
            const NumericVector& posterior_prob,
            const CharacterVector& col_names,
            double threshold,
            double r2_threshold,
            double coverage) {
    
    try {
        int p = mat.nrow();
        int L = mat.ncol();
        int m = cmfg_mat.nrow();
        std::vector<std::vector<int>> clusters;
        std::vector<double> mps_values;
        std::vector<std::vector<double>> r2_values;

        // Sort each column in descending order
        vector<vector<pair<double, int>>> sortedMat(L);
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < p; j++) {
                if (mat(j, i) >= threshold) {
                    sortedMat[i].emplace_back(mat(j, i), j); // Only include SNPs above threshold
                }
            }
            sort(sortedMat[i].begin(), sortedMat[i].end(), greater<pair<double, int>>());
        }

        // For each column in alpha matrix
        for (int col = 0; col < L; col++) {
            int starting_idx = 0;

            // If this is the last column and the null SNP is top, start from the second SNP
            if (col == L-1 && sortedMat[col][0].second == p - 1) {
                if (sortedMat[col].size() < 2) continue;
                starting_idx = 1;
            }

            int top_snp = sortedMat[col][starting_idx].second;
            if (top_snp < 0) continue;

            // Start with the top SNP
            std::vector<int> current_cluster = {top_snp};
            std::vector<double> current_r2s;

            // Check posterior probabilities
            double sum_prob_without = 0.0;
            std::vector<int> models_without_snps;
            
            for(int i = 0; i < m; i++) {
                if (cmfg_mat(i, top_snp) == 0) {
                    models_without_snps.push_back(i);
                    sum_prob_without += posterior_prob[i];
                }
            }
            double mps = 1.0 - sum_prob_without;

            bool cluster_done = false;
            if (coverage <= 1 && mps >= coverage) {
                cluster_done = true;
            } else {
                // Add more SNPs if needed
                for (size_t j = starting_idx + 1; j < sortedMat[col].size(); j++){
                    int next_snp = sortedMat[col][j].second;
                    if (next_snp == p - 1 || next_snp < 0) continue;

                    bool add_to_cluster = true;
                    double r2;

                    // Check R2 with all SNPs in current cluster
                    for (int cluster_snp : current_cluster) {
                        r2 = compute_r2(X, cluster_snp, next_snp);
                        if (r2 < r2_threshold) {
                            add_to_cluster = false;
                            break;
                        }
                        current_r2s.push_back(r2);
                    }

                    if (add_to_cluster) {
                        current_cluster.push_back(next_snp);

                        // Update MPS
                        for (size_t idx = 0; idx < models_without_snps.size();) {
                            int model_idx = models_without_snps[idx];
                            if (cmfg_mat(model_idx, next_snp) == 1) { // if the model contains the next SNP, remove it from models_without_snps
                                sum_prob_without -= posterior_prob[model_idx];
                                models_without_snps[idx] = models_without_snps.back();
                                models_without_snps.pop_back();
                            } else {
                                idx++; // go to the next model
                            }
                        }
                        mps = 1.0 - sum_prob_without;

                        if (coverage <= 1 && mps >= coverage) {
                            cluster_done = true;
                            break;
                        }
                    }
                }
            }

            if (coverage > 1 || cluster_done) {
                clusters.push_back(current_cluster);
                mps_values.push_back(mps);
                r2_values.push_back(current_r2s);
            }
        }

        // Convert results to R objects with safety checks
        int n_clusters = clusters.size();
        List r_clusters(n_clusters);
        NumericVector r_mps(n_clusters);
        IntegerVector cluster_sizes(n_clusters);
        List r2_stats(n_clusters);

        for(int i = 0; i < n_clusters; i++) {
            CharacterVector cluster_names(clusters[i].size());
            for(size_t j = 0; j < clusters[i].size(); j++) {
                cluster_names[j] = col_names[clusters[i][j]];
            }
            
            r_clusters[i] = cluster_names;
            r_mps[i] = mps_values[i];
            cluster_sizes[i] = clusters[i].size();

            vector<double>& r2s = r2_values[i];
            if (!r2s.empty()) {
                sort(r2s.begin(), r2s.end());
                r2_stats[i] = DataFrame::create(
                    Named("min_r2") = r2s.front(),
                    Named("mean_r2") = accumulate(r2s.begin(), r2s.end(), 0.0) / r2s.size(),
                    Named("median_r2") = median(r2s)
                );
            } else {
                r2_stats[i] = DataFrame::create(
                    Named("min_r2") = 1,
                    Named("mean_r2") = 1,
                    Named("median_r2") = 1
                );
            }
        }

        return List::create(
            Named("clusters") = r_clusters,
            Named("spip") = r_mps,
            Named("size") = cluster_sizes,
            Named("cluster_r2") = r2_stats,
            Named("r2_threshold") = r2_threshold,
            Named("coverage") = String(coverage > 1 ? "signal cluster" : std::to_string(static_cast<int>(coverage * 100)) + "% credible set")
        );
        
    } catch (std::exception& e) {
        Rcpp::stop("Error in get_sc: %s", e.what());
    } catch (...) {
        Rcpp::stop("Unknown error in get_sc");
    }
}