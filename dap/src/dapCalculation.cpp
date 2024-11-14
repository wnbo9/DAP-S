#include <Rcpp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <map>
#include <cmath>

using namespace Rcpp;
using namespace std;

class MLR {
public:
  void init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_);
  void set_effect_vec(const vector<double> &phi2_vec);
  double compute_log10_BF(const vector<int> &indicator);
  double log10_weighted_sum(const vector<double> &vec, const vector<double> &wts);

private:
  int n, p;
  double yty;
  gsl_matrix *GtG;
  gsl_matrix *Gty;
  vector<double> phi2_vec;
};

void MLR::init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_) {
  n = n_;
  p = GtG_->size1;
  yty = yty_;

  GtG = gsl_matrix_calloc(p, p);
  Gty = gsl_matrix_calloc(p, 1);
  gsl_matrix_memcpy(GtG, GtG_);
  gsl_matrix_memcpy(Gty, Gty_);
}

void MLR::set_effect_vec(const vector<double> &phi2_vec_) {
  phi2_vec = phi2_vec_;
}

double MLR::compute_log10_BF(const vector<int> &indicator) {
  vector<double> rstv;
  vector<double> wv;

  int ep = 0;
  int count = 0;
  map<int, int> imap;

  for (int i = 0; i < indicator.size(); i++) {
    if (indicator[i] == 1) {
      ep++;
      imap[i] = count++;
    }
  }

  if (ep == 0) {
    return 0;
  }

  gsl_matrix *XtX = gsl_matrix_calloc(ep, ep);
  gsl_matrix *Xty = gsl_matrix_calloc(ep, 1);

  for (int i = 0; i < p; i++) {
    if (indicator[i] == 0)
      continue;
    gsl_matrix_set(Xty, imap[i], 0, gsl_matrix_get(Gty, i, 0));
    for (int j = 0; j < p; j++) {
      if (indicator[j] == 1) {
        double val = gsl_matrix_get(GtG, i, j);
        gsl_matrix_set(XtX, imap[i], imap[j], val);
      }
    }
  }

  gsl_matrix *V = gsl_matrix_calloc(ep, ep);
  gsl_vector *S = gsl_vector_calloc(ep);
  gsl_vector *work = gsl_vector_calloc(ep);
  gsl_linalg_SV_decomp(XtX, V, S, work);

  for (int i = 0; i < phi2_vec.size(); i++) {
    double det = 1;
    double phi2 = phi2_vec[i];

    gsl_matrix *tt1 = gsl_matrix_calloc(ep, ep);

    for (int j = 0; j < ep; j++) {
      det = det * (1 + phi2 * gsl_vector_get(S, j));

      double v = gsl_vector_get(S, j) + 1.0 / phi2_vec[i];
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
    wv.push_back(1.0 / phi2_vec.size());

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

//' Compute log10 prior of model configuration
//' @param mcfg A vector of model configuration
//' @param pi_vec A vector of prior probabilities
//' @return Result Log 10 prior of model configuration
//' @export
// [[Rcpp::export]]
 double compute_log10_prior(IntegerVector mcfg, NumericVector pi_vec) {
   double lp = 0;
   int p = mcfg.size();

   for (int i = 0; i < p; i++) {
     if (mcfg[i] == 0) {
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
 NumericVector compute_log10_posterior(NumericMatrix X, NumericVector y, IntegerMatrix cmfg_matrix, NumericVector pi_vec, NumericVector phi2_vec) {
   int n = X.nrow();
   int p = X.ncol();
   int m = cmfg_matrix.nrow();

   // Compute GtG and Gty
   gsl_matrix *GtG = gsl_matrix_calloc(p, p);
   gsl_matrix *Gty = gsl_matrix_calloc(p, 1);
   double yty = 0;

   for (int i = 0; i < n; i++) {
     for (int j = 0; j < p; j++) {
       for (int k = 0; k < p; k++) {
         double val = gsl_matrix_get(GtG, j, k) + X(i, j) * X(i, k);
         gsl_matrix_set(GtG, j, k, val);
       }
       double val = gsl_matrix_get(Gty, j, 0) + X(i, j) * y[i];
       gsl_matrix_set(Gty, j, 0, val);
     }
     yty += y[i] * y[i];
   }

   // Initialize MLR
   MLR mlr;
   mlr.init(yty, GtG, Gty, n);
   mlr.set_effect_vec(as<vector<double>>(phi2_vec));

   // Compute log10 posterior scores for each model configuration
   NumericVector log10_posterior_scores(m);
   for (int i = 0; i < m; i++) {
     IntegerVector cmfg = cmfg_matrix(i, _);
     vector<int> indicator(cmfg.begin(), cmfg.end());

     double log10BF = mlr.compute_log10_BF(indicator);
     double log10_prior = compute_log10_prior(cmfg, pi_vec);

     log10_posterior_scores[i] = log10BF + log10_prior;
   }

   // Free allocated memory
   gsl_matrix_free(GtG);
   gsl_matrix_free(Gty);

   return log10_posterior_scores;
 }
