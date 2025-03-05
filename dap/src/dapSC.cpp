#include <Rcpp.h>
#include <gsl/gsl_statistics.h>
#include <vector>
#include <set>
#include <algorithm>

using namespace Rcpp;
using namespace std;


double compute_r2(const NumericMatrix& X, int i, int j) {
    double r2;
    int n = X.nrow();
    const double* gi = &X(0, i);
    const double* gj = &X(0, j);
    
    r2 = pow(gsl_stats_correlation(gi, 1, gj, 1, n), 2);
    return r2;
}

double compute_r2_ss(const NumericMatrix& XtX, int i, int j) {
    // For diagonal elements, r2 is 1
    if (i == j) return 1.0;
    
    // Extract correlation from XtX matrix and square it
    double corr = XtX(i, j) / sqrt(XtX(i, i) * XtX(j, j));
    double r2 = corr * corr;
    
    return r2;
}

double median(vector<double>& v) {
    size_t n = v.size();
    if(n % 2 == 0) {
        return (v[n/2 - 1] + v[n/2]) / 2;
    }
    return v[n/2];
}

//' Get signal clusters or credible sets at given coverage level
//' @param X NumericMatrix containing the raw data
//' @param effect_pip NumericMatrix containing the effect posterior inclusion probabilities
//' @param snp_names CharacterVector of SNP names
//' @param r2_threshold Double specifying the R-squared threshold for correlation
//' @param coverage Double specifying the coverage threshold
//' @return List containing:
//'   \item{clusters}{List of character vectors containing cluster memberships}
//'   \item{spip}{Numeric vector of signal posterior inclusion probabilities}
//'   \item{sizes}{Integer vector of cluster sizes}
//'   \item{r2_threshold}{Double containing the R-squared threshold used}
//'   \item{coverage}{Double containing the coverage threshold used}
List get_sc(const NumericMatrix& X,
            const NumericMatrix& effect_pip,
            const CharacterVector& snp_names,
            double r2_threshold,
            double coverage) {
    
    try {
        int p = X.ncol(); // 5000
        int L = effect_pip.ncol(); // 3

        vector<vector<int>> clusters;
        vector<double> cluster_pips;
        vector<vector<double>> r2_values;
        vector<int> sc_index;

        // For each column
        for (int col = 0; col < L; col++) {
            vector<pair<double, int>> snp_probs;
            for (int snp = 0; snp < p; snp++) {
                double pip = effect_pip(snp, col);
                if (pip > 1e-6) {
                    snp_probs.emplace_back(pip, snp);
                }
            }

            // Sort SNPs by their PIPs in descending order
            sort(snp_probs.begin(), snp_probs.end(), greater<pair<double, int>>());
            if (snp_probs.empty()) continue;

            // Start with the top SNP
            int top_snp = snp_probs[0].second;
            vector<int> current_cluster = {top_snp};
            vector<double> current_r2s;
            double current_sum_pip = snp_probs[0].first;

            bool cluster_done = false;
            if (coverage <= 1 && current_sum_pip >= coverage) {
                cluster_done = true;
            } else {
                // Process remaining SNPs in the rank order
                for (size_t j = 1; j < snp_probs.size(); j++){
                    int next_snp = snp_probs[j].second;
                    bool add_to_cluster = true;

                    // Check R2 with all SNPs in current cluster
                    for (int cluster_snp : current_cluster) {
                        double r2 = compute_r2(X, cluster_snp, next_snp);
                        if (r2 < r2_threshold) {
                            add_to_cluster = false;
                            break;
                        }
                        current_r2s.push_back(r2);
                    }

                    if (add_to_cluster) {
                        current_cluster.push_back(next_snp);
                        current_sum_pip += snp_probs[j].first;

                        if (coverage <= 1 && current_sum_pip >= coverage) {
                            cluster_done = true;
                            break;
                        }
                    }
                }
            }

            if (coverage > 1 || cluster_done) {
                clusters.push_back(current_cluster);
                cluster_pips.push_back(max(0.0, min(1.0, round(current_sum_pip * 1000.0) / 1000.0)));
                r2_values.push_back(current_r2s);
                sc_index.push_back(col+1);
            }
        }

        // Convert results to R objects with safety checks
        int n_clusters = clusters.size();
        List r_clusters(n_clusters);
        NumericVector r_mps(n_clusters);
        IntegerVector cluster_sizes(n_clusters);

        NumericVector min_r2(n_clusters);
        NumericVector mean_r2(n_clusters);
        NumericVector median_r2(n_clusters);

        for(int i = 0; i < n_clusters; i++) {
            CharacterVector cluster_names(clusters[i].size());

            for(size_t j = 0; j < clusters[i].size(); j++) {
                cluster_names[j] = snp_names[clusters[i][j]];
            }
            
            r_clusters[i] = cluster_names;
            r_mps[i] = cluster_pips[i];
            cluster_sizes[i] = clusters[i].size();

            vector<double>& r2s = r2_values[i];
            if (!r2s.empty()) {
                sort(r2s.begin(), r2s.end());
                min_r2[i] = r2s.front();
                mean_r2[i] = accumulate(r2s.begin(), r2s.end(), 0.0) / r2s.size();
                median_r2[i] = median(r2s);
            } else {
                min_r2[i] = 1.0;
                mean_r2[i] = 1.0;
                median_r2[i] = 1.0;
            }
        }

        CharacterVector cluster_names_vec(n_clusters);
        for(int i = 0; i < n_clusters; i++) {
            cluster_names_vec[i] = "C" + to_string(sc_index[i]);
        }
        DataFrame r2_stats = DataFrame::create(
            Named("min_r2") = min_r2,
            Named("mean_r2") = mean_r2,
            Named("median_r2") = median_r2
        );
        r2_stats.attr("row.names") = cluster_names_vec;

        vector<vector<int>> snp_index = clusters;
        for(auto& cluster : snp_index) {
            transform(cluster.begin(), cluster.end(), cluster.begin(), [](int x) { return x + 1; });
        }

        List named_clusters;
        List named_snp_index;
        for(int i = 0; i < n_clusters; i++) {
            String cluster_name = as<std::string>(cluster_names_vec[i]);
            named_clusters[cluster_name] = r_clusters[i];
            named_snp_index[cluster_name] = snp_index[i];
        }

        
        return List::create(
            Named("cluster") = named_clusters,
            Named("snp_index") = named_snp_index,
            Named("cpip") = r_mps,
            Named("size") = cluster_sizes,
            Named("cluster_r2") = r2_stats,
            Named("sc_index") = sc_index,
            Named("r2_threshold") = r2_threshold,
            Named("requested_coverage") = String(coverage > 1 ? "signal cluster" : to_string(static_cast<int>(coverage * 100)) + "% credible set")
        );
        
    } catch (std::exception& e) {
        Rcpp::stop("Error in get_sc: %s", e.what());
    } catch (...) {
        Rcpp::stop("Unknown error in get_sc");
    }
}




//' Get signal clusters or credible sets at given coverage level
//' @param XtX NumericMatrix containing the raw data
//' @param effect_pip NumericMatrix containing the effect posterior inclusion probabilities
//' @param snp_names CharacterVector of SNP names
//' @param r2_threshold Double specifying the R-squared threshold for correlation
//' @param coverage Double specifying the coverage threshold
//' @return List containing:
//'   \item{clusters}{List of character vectors containing cluster memberships}
//'   \item{spip}{Numeric vector of signal posterior inclusion probabilities}
//'   \item{sizes}{Integer vector of cluster sizes}
//'   \item{r2_threshold}{Double containing the R-squared threshold used}
//'   \item{coverage}{Double containing the coverage threshold used}
List get_sc_ss(const NumericMatrix& XtX,
    const NumericMatrix& effect_pip,
    const CharacterVector& snp_names,
    double r2_threshold,
    double coverage) {

try {
int p = XtX.ncol(); // 5000
int L = effect_pip.ncol(); // 3

vector<vector<int>> clusters;
vector<double> cluster_pips;
vector<vector<double>> r2_values;
vector<int> sc_index;

// For each column
for (int col = 0; col < L; col++) {
    vector<pair<double, int>> snp_probs;
    for (int snp = 0; snp < p; snp++) {
        double pip = effect_pip(snp, col);
        if (pip > 1e-6) {
            snp_probs.emplace_back(pip, snp);
        }
    }

    // Sort SNPs by their PIPs in descending order
    sort(snp_probs.begin(), snp_probs.end(), greater<pair<double, int>>());
    if (snp_probs.empty()) continue;

    // Start with the top SNP
    int top_snp = snp_probs[0].second;
    vector<int> current_cluster = {top_snp};
    vector<double> current_r2s;
    double current_sum_pip = snp_probs[0].first;

    bool cluster_done = false;
    if (coverage <= 1 && current_sum_pip >= coverage) {
        cluster_done = true;
    } else {
        // Process remaining SNPs in the rank order
        for (size_t j = 1; j < snp_probs.size(); j++){
            int next_snp = snp_probs[j].second;
            bool add_to_cluster = true;

            // Check R2 with all SNPs in current cluster
            for (int cluster_snp : current_cluster) {
                double r2 = compute_r2_ss(XtX, cluster_snp, next_snp);
                if (r2 < r2_threshold) {
                    add_to_cluster = false;
                    break;
                }
                current_r2s.push_back(r2);
            }

            if (add_to_cluster) {
                current_cluster.push_back(next_snp);
                current_sum_pip += snp_probs[j].first;

                if (coverage <= 1 && current_sum_pip >= coverage) {
                    cluster_done = true;
                    break;
                }
            }
        }
    }

    if (coverage > 1 || cluster_done) {
        clusters.push_back(current_cluster);
        cluster_pips.push_back(max(0.0, min(1.0, round(current_sum_pip * 1000.0) / 1000.0)));
        r2_values.push_back(current_r2s);
        sc_index.push_back(col+1);
    }
}

// Convert results to R objects with safety checks
int n_clusters = clusters.size();
List r_clusters(n_clusters);
NumericVector r_mps(n_clusters);
IntegerVector cluster_sizes(n_clusters);

NumericVector min_r2(n_clusters);
NumericVector mean_r2(n_clusters);
NumericVector median_r2(n_clusters);

for(int i = 0; i < n_clusters; i++) {
    CharacterVector cluster_names(clusters[i].size());

    for(size_t j = 0; j < clusters[i].size(); j++) {
        cluster_names[j] = snp_names[clusters[i][j]];
    }
    
    r_clusters[i] = cluster_names;
    r_mps[i] = cluster_pips[i];
    cluster_sizes[i] = clusters[i].size();

    vector<double>& r2s = r2_values[i];
    if (!r2s.empty()) {
        sort(r2s.begin(), r2s.end());
        min_r2[i] = r2s.front();
        mean_r2[i] = accumulate(r2s.begin(), r2s.end(), 0.0) / r2s.size();
        median_r2[i] = median(r2s);
    } else {
        min_r2[i] = 1.0;
        mean_r2[i] = 1.0;
        median_r2[i] = 1.0;
    }
}

CharacterVector cluster_names_vec(n_clusters);
for(int i = 0; i < n_clusters; i++) {
    cluster_names_vec[i] = "C" + to_string(sc_index[i]);
}
DataFrame r2_stats = DataFrame::create(
    Named("min_r2") = min_r2,
    Named("mean_r2") = mean_r2,
    Named("median_r2") = median_r2
);
r2_stats.attr("row.names") = cluster_names_vec;

vector<vector<int>> snp_index = clusters;
for(auto& cluster : snp_index) {
    transform(cluster.begin(), cluster.end(), cluster.begin(), [](int x) { return x + 1; });
}

List named_clusters;
List named_snp_index;
for(int i = 0; i < n_clusters; i++) {
    String cluster_name = as<std::string>(cluster_names_vec[i]);
    named_clusters[cluster_name] = r_clusters[i];
    named_snp_index[cluster_name] = snp_index[i];
}


return List::create(
    Named("cluster") = named_clusters,
    Named("snp_index") = named_snp_index,
    Named("cpip") = r_mps,
    Named("size") = cluster_sizes,
    Named("cluster_r2") = r2_stats,
    Named("sc_index") = sc_index,
    Named("r2_threshold") = r2_threshold,
    Named("requested_coverage") = String(coverage > 1 ? "signal cluster" : to_string(static_cast<int>(coverage * 100)) + "% credible set")
);

} catch (std::exception& e) {
Rcpp::stop("Error in get_sc: %s", e.what());
} catch (...) {
Rcpp::stop("Unknown error in get_sc");
}
}