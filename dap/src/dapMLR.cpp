#include <Rcpp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

using namespace Rcpp;
using namespace std;

class MLR {
private:
    int n, p;
    double yty;
    gsl_matrix *GtG;
    gsl_matrix *Gty;
    NumericMatrix phi2_mat;

public:
    void init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_, const NumericMatrix &phi2_mat_);
    double compute_log10_BF(const vector<int> &indicator);
    double log10_weighted_sum(const vector<double> &vec, const vector<double> &wts);

    ~MLR() {
        if (GtG != 0)
            gsl_matrix_free(GtG);
        if (Gty != 0)
            gsl_matrix_free(Gty);
    }
};

void MLR::init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_, const NumericMatrix &phi2_mat_) {
  n = n_;
  p = GtG_->size1;
  yty = yty_;
  phi2_mat = phi2_mat_;

  GtG = gsl_matrix_calloc(p, p);
  Gty = gsl_matrix_calloc(p, 1);
  gsl_matrix_memcpy(GtG, GtG_);
  gsl_matrix_memcpy(Gty, Gty_);
}

double MLR::compute_log10_BF(const vector<int> &indicator) {
  vector<double> rstv;
  vector<double> wv;

  int ep = 0;
  int count = 0;
  map<int, int> imap;

  for (int i = 0; i < indicator.size(); i++) {
      if (indicator[i] != p) {
        ep++;
        imap[i] = count++;
      }
  }

  if (ep == 0) {
      return 0;
  }

  gsl_matrix *XtX = gsl_matrix_calloc(ep, ep);
  gsl_matrix *Xty = gsl_matrix_calloc(ep, 1);

  for (int i = 0; i < indicator.size(); i++) {
      if (indicator[i] == p)
          continue;

      gsl_matrix_set(Xty, imap[i], 0, gsl_matrix_get(Gty, indicator[i], 0));
      
      for (int j = 0; j < indicator.size(); j++) {
        if (indicator[j] != p) {
          double val = gsl_matrix_get(GtG, indicator[i], indicator[j]);
          gsl_matrix_set(XtX, imap[i], imap[j], val);
        }
      }
  }

  gsl_matrix *V = gsl_matrix_calloc(ep, ep);
  gsl_vector *S = gsl_vector_calloc(ep);
  gsl_vector *work = gsl_vector_calloc(ep);
  gsl_linalg_SV_decomp(XtX, V, S, work);

  for (int i = 0; i < phi2_mat.nrow(); i++) {
      double det = 1;

      gsl_matrix *tt1 = gsl_matrix_calloc(ep, ep);

      for (int j = 0; j < ep; j++) {
        det = det * (1 + phi2_mat(i, imap[j]) * gsl_vector_get(S, j));

        double v = gsl_vector_get(S, j) + 1.0 / phi2_mat(i, imap[j]);
        if (v > 1e-8) {
          gsl_matrix_set(tt1, j, j, 1.0 / v);
        }
      }

      gsl_matrix *tt2 = gsl_matrix_calloc(ep, ep);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, tt1, 0, tt2);

      gsl_matrix *IWV = gsl_matrix_calloc(ep, ep);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tt2, V, 0, IWV);

      gsl_matrix *tt3 = gsl_matrix_calloc(ep, 1);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, IWV, Xty, 0, tt3);

      gsl_matrix *tt4 = gsl_matrix_calloc(1, 1);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Xty, tt3, 0, tt4);

      double b_rss = gsl_matrix_get(tt4, 0, 0);
      double log10BF = -0.5 * log10(fabs(det)) - 0.5 * n * log10(1 - b_rss / yty);

      rstv.push_back(log10BF);
      wv.push_back(1.0 / phi2_mat.nrow());

      gsl_matrix_free(tt1);
      gsl_matrix_free(tt2);
      gsl_matrix_free(IWV);
      gsl_matrix_free(tt3);
      gsl_matrix_free(tt4);
  }

  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);

  gsl_matrix_free(XtX);
  gsl_matrix_free(Xty);

  double rst = log10_weighted_sum(rstv, wv);
  
  return rst;
}

  
double MLR::log10_weighted_sum(const vector<double> &vec, const vector<double> &wts) {
    double max = vec[0];
    for (size_t i = 0; i < vec.size(); i++) {
        if (vec[i] > max)
            max = vec[i];
    }
    double sum = 0;
    for (size_t i = 0; i < vec.size(); i++) {
        sum += wts[i] * pow(10, (vec[i] - max));
    }
    return (max + log10(sum));
}

double compute_log10_prior(const std::vector<int> &mcfg, NumericVector pi_vec) {
    double lp = 0;
    int p = pi_vec.length();
    int L = mcfg.size();

    std::vector<int> model_indices(p+1, 0);
    for (int i = 0; i < L; i++) {
        model_indices[mcfg[i]] = 1;
    }
    for (int i = 0; i < p; i++) {
        if (model_indices[i] == 0) {
            lp += log(1 - pi_vec[i]);
        } else {
            lp += log(pi_vec[i]);
        }
    }
    return lp / log(10);
}

//' Compute log10 posterior scores for multiple model configurations
//' @param X An n*p matrix of genotype data
//' @param y An n-vector of phenotype data
//' @param cmfg_matrix An m*p matrix of model configurations
//' @param pi_vec A vector of prior probabilities
//' @param phi2_vec A vector of phi2 values
//' @return An m-vector of log10 posterior scores
//' @export
// [[Rcpp::export]]
List compute_log10_posterior(const NumericMatrix& X, const NumericVector& y, 
                             const std::vector<std::vector<int>>& cmfg_matrix,
                             const NumericVector& pi_vec, const NumericMatrix& phi2_mat) {
  int n = X.nrow();
  int p = X.ncol();
  int m = cmfg_matrix.size();

  // Rest of initialization code...
  gsl_matrix *Y = gsl_matrix_calloc(n, 1);
  double yty = 0;
  for (int i = 0; i < n; i++) {
    double val = y[i];
    yty += val * val;
    gsl_matrix_set(Y, i, 0, val);
  }

  gsl_matrix *G = gsl_matrix_calloc(n, p);
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      gsl_matrix_set(G, i, j, X(i, j));
    }
  }

  gsl_matrix *GtG = gsl_matrix_calloc(p, p);
  gsl_matrix *Gty = gsl_matrix_calloc(p, 1);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, G, G, 0.0, GtG);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, G, Y, 0.0, Gty);

  MLR mlr;
  mlr.init(yty, GtG, Gty, n, phi2_mat);
  
  NumericVector log10_BF(m);
  NumericVector log10_prior(m);
  NumericVector log10_posterior_score(m);

  for (int i = 0; i < m; i++) {
    // Extract the row
    vector<int> indicator = cmfg_matrix[i];

    log10_BF[i] = mlr.compute_log10_BF(indicator);
    log10_prior[i] = compute_log10_prior(indicator, pi_vec);
    log10_posterior_score[i] = log10_BF[i] + log10_prior[i];
  }

  // Free allocated memory
  gsl_matrix_free(G);
  gsl_matrix_free(Y);
  gsl_matrix_free(GtG);
  gsl_matrix_free(Gty);

  return List::create(
    Named("log10_BF") = log10_BF,
    Named("log10_prior") = log10_prior,
    Named("log10_posterior_score") = log10_posterior_score
  );
}