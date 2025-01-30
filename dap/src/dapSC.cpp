#include <Rcpp.h>
#include <gsl/gsl_statistics.h>
#include <vector>
#include <set>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(gsl)]]

double compute_r2(const NumericMatrix& X, int i, int j) {
    double r2;
    int n = X.nrow();
    const double* gi = &X(0, i);
    const double* gj = &X(0, j);
    
    r2 = pow(gsl_stats_correlation(gi, 1, gj, 1, n), 2);
    return r2;
}

double median(vector<double>& v) {
    size_t n = v.size();
    if(n % 2 == 0) {
        return (v[n/2 - 1] + v[n/2]) / 2;
    }
    return v[n/2];
}


double calculate_mps(const std::vector<int>& cluster, 
                    const NumericMatrix& combo,
                    const NumericMatrix& single,
                    const NumericVector& posterior_prob) {
    double mps = 0.0;
    int m1 = combo.nrow();
    int m2 = single.nrow();
    int L = combo.ncol();
    
    // Check combo models
    for (int i = 0; i < m1; i++) {
        bool includes_cluster_snp = false;
        for (int col = 0; col < L; col++) {
            int model_snp = combo(i, col);
            for (int cluster_snp : cluster) {
                if (model_snp == cluster_snp) {
                    includes_cluster_snp = true;
                    break;
                }
            }
            if (includes_cluster_snp) break;
        }
        if (includes_cluster_snp) {
            mps += posterior_prob[i];
        }
    }
    
    // Check single SNP models
    for (int i = 0; i < m2; i++) {
        int model_snp = single(i, 0);
        for (int cluster_snp : cluster) {
            if (model_snp == cluster_snp) {
                mps += posterior_prob[m1 + i];
                break;
            }
        }
    }
    
    return mps;
}

//' Get signal clusters or credible sets at given coverage level
//' 
//' This function identifies clusters of correlated signals based on R-squared values
//' and model posterior probabilities.
//' 
//' @param X NumericMatrix containing the raw data
//' @param combo NumericMatrix containing the model configurations
//' @param single NumericMatrix containing the single SNP models
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
            const NumericMatrix& combo,
            const NumericMatrix& single,
            const NumericVector& posterior_prob,
            const CharacterVector& col_names,
            double threshold,
            double r2_threshold,
            double coverage) {
    
    try {
        int p = X.ncol(); // 5000
        int L = combo.ncol(); // 3
        int m1 = combo.nrow();

        // Create PIP matrix for each column (L+1 columns: L columns)
        NumericMatrix pip_matrix(p, L);
        std::vector<std::vector<int>> clusters;
        std::vector<double> mps_values;
        std::vector<std::vector<double>> r2_values;
        std::vector<std::vector<double>> snp_pips;

        // For each column in alpha matrix
        for (int col = 0; col < L; col++) {
            
            vector<pair<double, int>> snp_probs;
            for (int snp = 0; snp < p; snp++) {
                double snp_prob = 0.0;
                // Check combo models for this column
                for (int model = 0; model < m1; model++) {
                    if (combo(model, col) == snp) {
                        snp_prob += posterior_prob[model];
                    }
                }
                pip_matrix(snp, col) = snp_prob;
                if (snp_prob>0) snp_probs.emplace_back(snp_prob, snp);
            }

            // Sort SNPs by their probabilities
            sort(snp_probs.begin(), snp_probs.end(), greater<pair<double, int>>());
            if (snp_probs.empty()) continue;

            // Start with the top SNP
            int top_snp = snp_probs[0].second;
            std::vector<int> current_cluster = {top_snp};
            std::vector<double> current_r2s;

            // Calculate initial MPS
            double mps = calculate_mps(current_cluster, combo, single, posterior_prob);

            bool cluster_done = false;
            if (coverage <= 1 && mps >= coverage) {
                cluster_done = true;
            } else {
                // Process remaining SNPs in the rank order
                for (size_t j = 1; j < snp_probs.size(); j++){
                    int next_snp = snp_probs[j].second;
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
                        mps = calculate_mps(current_cluster, combo, single, posterior_prob);

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

        std::vector<std::vector<int>> snp_index = clusters;
        for(auto& cluster : snp_index) {
            std::transform(cluster.begin(), cluster.end(), cluster.begin(), [](int x) { return x + 1; });
        }

        return List::create(
            Named("clusters") = r_clusters,
            Named("snp_index") = snp_index,
            Named("spip") = r_mps,
            Named("size") = cluster_sizes,
            Named("cluster_r2") = r2_stats,
            Named("pip_matrix") = pip_matrix,
            Named("r2_threshold") = r2_threshold,
            Named("coverage") = String(coverage > 1 ? "signal cluster" : std::to_string(static_cast<int>(coverage * 100)) + "% credible set")
        );
        
    } catch (std::exception& e) {
        Rcpp::stop("Error in get_sc: %s", e.what());
    } catch (...) {
        Rcpp::stop("Unknown error in get_sc");
    }
}