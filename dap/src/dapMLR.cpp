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
    vector<double> compute_log10_BF(const vector<int> &indicator, bool twas_weight);
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

vector<double> MLR::compute_log10_BF(const vector<int> &indicator, bool twas_weight = true) {
  vector<double> rstv;
  vector<double> wv;
  vector<double> result(p+1, 0.0);

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
    return result;
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

  // Calculate log10 BF
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

  // Calculate regression weights
  if (twas_weight) {
    gsl_matrix *t1 = gsl_matrix_calloc(ep, ep);
    for (int j = 0; j < ep; j++) {
        double v = gsl_vector_get(S, j);
        if (v > 1e-8) {
          gsl_matrix_set(t1, j, j, 1.0 / v);
        }
    }
    gsl_matrix *t2 = gsl_matrix_calloc(ep,ep);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, t1, 0, t2);

    gsl_matrix *XtX_inv = gsl_matrix_calloc(ep,ep);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, t2, V, 0, XtX_inv);

    // (X'X)^{-1)X'y
    gsl_matrix *t3 = gsl_matrix_calloc(ep,1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, XtX_inv, Xty, 0, t3);

    for (int j = 0; j < indicator.size(); j++) {
        if (indicator[j] == p)
            continue;
        int pos = indicator[j];
        result[pos + 1] = gsl_matrix_get(t3, imap[j], 0);
    }

    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
    gsl_matrix_free(XtX_inv);
    gsl_matrix_free(t3);
  }


  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);

  gsl_matrix_free(XtX);
  gsl_matrix_free(Xty);

  result[0] = log10_weighted_sum(rstv, wv);
  
  return result;
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


// *** NEW: Class for infinitesimal model that mirrors MLR structure ***
class MLR_infinitesimal {
private:
    int p;
    gsl_matrix *XtOmegaX;  // *** NEW: Pre-computed for full X (like GtG in MLR) ***
    gsl_vector *XtOmegay;  // *** NEW: Pre-computed for full X (like Gty in MLR) ***
    NumericMatrix phi2_mat;

public:
    void init(const gsl_matrix *XtOmegaX_, const gsl_vector *XtOmegay_, const NumericMatrix &phi2_mat_);
    vector<double> compute_log10_BF(const vector<int> &indicator, bool twas_weight);
    double log10_weighted_sum(const vector<double> &vec, const vector<double> &wts);  // *** ADDED: Missing function ***
    
    ~MLR_infinitesimal() {
        if (XtOmegaX != 0)
            gsl_matrix_free(XtOmegaX);
        if (XtOmegay != 0)
            gsl_vector_free(XtOmegay);
    }
};

void MLR_infinitesimal::init(const gsl_matrix *XtOmegaX_, const gsl_vector *XtOmegay_, const NumericMatrix &phi2_mat_) {
    p = XtOmegaX_->size1;
    phi2_mat = phi2_mat_;
    
    XtOmegaX = gsl_matrix_calloc(p, p);
    XtOmegay = gsl_vector_calloc(p);
    gsl_matrix_memcpy(XtOmegaX, XtOmegaX_);
    gsl_vector_memcpy(XtOmegay, XtOmegay_);
}

vector<double> MLR_infinitesimal::compute_log10_BF(const vector<int> &indicator, bool twas_weight = false) {
    vector<double> rstv;
    vector<double> wv;
    vector<double> result(p+1, 0.0);
    
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
        return result;
    }
    
    // *** CHANGED: Extract submatrix and subvector for the model (similar to MLR) ***
    gsl_matrix *model_XtOmegaX = gsl_matrix_calloc(ep, ep);
    gsl_vector *model_XtOmegay = gsl_vector_calloc(ep);
    
    // *** CHANGED: Extract from pre-computed full matrices instead of computing from scratch ***
    for (int i = 0; i < indicator.size(); i++) {
        if (indicator[i] == p)
            continue;
            
        gsl_vector_set(model_XtOmegay, imap[i], gsl_vector_get(XtOmegay, indicator[i]));
        
        for (int j = 0; j < indicator.size(); j++) {
            if (indicator[j] != p) {
                double val = gsl_matrix_get(XtOmegaX, indicator[i], indicator[j]); // *** Extract from pre-computed ***
                gsl_matrix_set(model_XtOmegaX, imap[i], imap[j], val);
            }
        }
    }
    
    // Compute SVD once for this model
    gsl_matrix *U = gsl_matrix_calloc(ep, ep);
    gsl_matrix_memcpy(U, model_XtOmegaX);
    gsl_matrix *V_svd = gsl_matrix_calloc(ep, ep);
    gsl_vector *S = gsl_vector_calloc(ep);
    gsl_vector *work = gsl_vector_calloc(ep);
    gsl_linalg_SV_decomp(U, V_svd, S, work);
    
    // For each phi2 value
    for (int t = 0; t < phi2_mat.nrow(); t++) {
        double s2 = phi2_mat(t, 0);
        double det = 1.0;
        
        // Build diagonal matrix for inverse computation (like MLR builds tt1)
        gsl_matrix* tt1 = gsl_matrix_calloc(ep, ep);
        for (int j = 0; j < ep; j++) {
            double s_j = gsl_vector_get(S, j);
            det *= (1.0 + s2 * s_j);
            
            double v = s_j + 1.0 / s2;
            if (v > 1e-8) {
                gsl_matrix_set(tt1, j, j, 1.0 / v);
            }
        }
        
        // Build the inverse matrix (like MLR builds IWV)
        gsl_matrix* omega_inv = gsl_matrix_calloc(ep, ep);
        gsl_matrix* temp2 = gsl_matrix_calloc(ep, ep);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_svd, tt1, 0.0, temp2);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp2, V_svd, 0.0, omega_inv);
        
        // Calculate quadratic form (like MLR computes with tt3 and tt4)
        gsl_vector* temp_vec = gsl_vector_calloc(ep);
        gsl_blas_dgemv(CblasNoTrans, 1.0, omega_inv, model_XtOmegay, 0.0, temp_vec);
        double quad_form = 0.0;
        gsl_blas_ddot(model_XtOmegay, temp_vec, &quad_form);
        
        // Calculate log Bayes Factor
        double log_bf = -0.5 * log(det) + 0.5 * quad_form;
        double log10_bf = log_bf / log(10.0);
        
        rstv.push_back(log10_bf);
        wv.push_back(1.0 / phi2_mat.nrow());
        
        // Free matrices allocated in loop (like MLR)
        gsl_matrix_free(tt1);
        gsl_matrix_free(temp2);
        gsl_matrix_free(omega_inv);
        gsl_vector_free(temp_vec);
    }
    
    // *** CHANGED: Use log10_weighted_sum like MLR does ***
    result[0] = log10_weighted_sum(rstv, wv);
    
    if (twas_weight) {
        // Compute beta_hat = (X^T Ω X)^{-1} X^T Ω y using SVD (like MLR)
        gsl_matrix *t1 = gsl_matrix_calloc(ep, ep);
        for (int j = 0; j < ep; j++) {
            double s_j = gsl_vector_get(S, j);
            if (s_j > 1e-8) {
                gsl_matrix_set(t1, j, j, 1.0 / s_j);
            }
        }
        
        gsl_matrix *t2 = gsl_matrix_calloc(ep, ep);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_svd, t1, 0.0, t2);
        
        gsl_matrix *XtOmegaX_inv = gsl_matrix_calloc(ep, ep);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, t2, V_svd, 0.0, XtOmegaX_inv);
        
        // Compute beta = (X^T Ω X)^{-1} X^T Ω y
        gsl_vector *beta_hat = gsl_vector_calloc(ep);
        gsl_blas_dgemv(CblasNoTrans, 1.0, XtOmegaX_inv, model_XtOmegay, 0.0, beta_hat);
        
        // Store the results in the appropriate slots
        for (int j = 0; j < indicator.size(); j++) {
            if (indicator[j] == p)
                continue;
            
            int pos = indicator[j];
            double weight = gsl_vector_get(beta_hat, imap[j]);
            result[pos + 1] = weight;
        }
        
        gsl_matrix_free(t1);
        gsl_matrix_free(t2);
        gsl_matrix_free(XtOmegaX_inv);
        gsl_vector_free(beta_hat);
    }
    
    // Free allocated memory
    gsl_matrix_free(model_XtOmegaX);
    gsl_vector_free(model_XtOmegay);
    gsl_matrix_free(U);
    gsl_matrix_free(V_svd);
    gsl_vector_free(S);
    gsl_vector_free(work);
    
    return result;
}

double MLR_infinitesimal::log10_weighted_sum(const vector<double> &vec, const vector<double> &wts) {
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
//' @param cmfg_matrix An m*p matrix of model configurations
//' @param pi_vec A vector of prior probabilities
//' @param phi2_mat A matrix of phi2 values
//' @param ss Flag indicating use of summary statistics (0 = use X,y; 1 = use summary statistics)
//' @param X_input An n*p matrix of genotype data (required when ss=0)
//' @param y_input An n-vector of phenotype data (required when ss=0)
//' @param XtX_input A p*p matrix of precomputed X'X (required when ss=1)
//' @param Xty_input A p*1 vector of precomputed X'y (required when ss=1)
//' @param yty_input A scalar representing y'y (required when ss=1)
//' @param n_input Sample size (required when ss=1)
//' @param twas_weight A boolean indicating whether to compute TWAS weights
//' @return An m-vector of log10 posterior scores
//' @export
// [[Rcpp::export]]
List compute_log10_posterior(
  const std::vector<std::vector<int>>& cmfg_matrix, 
  const NumericVector& pi_vec, 
  const NumericMatrix& phi2_mat, 
  int ss = 0, 
  SEXP X_input = R_NilValue, 
  SEXP y_input = R_NilValue,
  SEXP XtX_input = R_NilValue, 
  SEXP Xty_input = R_NilValue, 
  SEXP yty_input = R_NilValue, 
  SEXP n_input = R_NilValue,
  SEXP V_input = R_NilValue,
  SEXP Dsq_input = R_NilValue,
  SEXP var_input = R_NilValue,
  SEXP XtOmegay_input = R_NilValue,
  bool twas_weight = false) {

  int n = 0, p = 0, m = cmfg_matrix.size();
  double yty = 0.0;
  gsl_matrix *GtG = NULL, *Gty = NULL;
  gsl_matrix *XtOmegaX_full = NULL;  // *** NEW: Pre-computed matrix for infinitesimal model ***
  gsl_vector *XtOmegay_full = NULL;  // *** NEW: Pre-computed vector for infinitesimal model ***
  NumericMatrix V;
  NumericVector Dsq, var, XtOmegay;

  if (ss == 0) {
    // Raw data mode (X and y)    
    NumericMatrix X = as<NumericMatrix>(X_input);
    NumericVector y = as<NumericVector>(y_input);
    
    n = X.nrow();
    p = X.ncol();
    
    gsl_matrix *Y = gsl_matrix_calloc(n, 1);
    yty = 0;
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

    GtG = gsl_matrix_calloc(p, p);
    Gty = gsl_matrix_calloc(p, 1);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, G, G, 0.0, GtG);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, G, Y, 0.0, Gty);

    gsl_matrix_free(G);
    gsl_matrix_free(Y);

  } else if (ss == 1) {
    // Summary statistics mode  
    NumericMatrix XtX = as<NumericMatrix>(XtX_input);
    NumericVector Xty = as<NumericVector>(Xty_input);
    yty = as<double>(yty_input);
    n = as<int>(n_input);
  
    p = XtX.ncol();
  
    // Convert summary statistics to GSL format
    GtG = gsl_matrix_calloc(p, p);
    Gty = gsl_matrix_calloc(p, 1);
  
    for (int i = 0; i < p; i++) {
      gsl_matrix_set(Gty, i, 0, Xty[i]);
      for (int j = 0; j < p; j++) {
        gsl_matrix_set(GtG, i, j, XtX(i, j));
      }
    }
  } else if (ss == 2) {
    // *** CHANGED: Infinitesimal model mode - now pre-compute XtOmegaX for full X ***
    V = as<NumericMatrix>(V_input);
    Dsq = as<NumericVector>(Dsq_input);
    var = as<NumericVector>(var_input);
    XtOmegay = as<NumericVector>(XtOmegay_input);
    p = Dsq.length() - 1;
    
    // *** NEW: Pre-compute X^T Ω X for the full X matrix (like GtG in MLR) ***
    XtOmegaX_full = gsl_matrix_calloc(p, p);
    XtOmegay_full = gsl_vector_calloc(p);
    
    // Copy XtOmegay to GSL format
    for (int i = 0; i < p; i++) {
        gsl_vector_set(XtOmegay_full, i, XtOmegay[i]);
    }
    
    // *** NEW: Compute X^T Ω X = V * (Dsq/var) * V^T for full matrix ***
    gsl_matrix* DsqVarInv = gsl_matrix_calloc(p, p);
    for (int i = 0; i < p; i++) {
        double dsq = Dsq[i];
        double v = var[i];
        double ratio = (v > 1e-8) ? dsq / v : dsq * 1e8;
        gsl_matrix_set(DsqVarInv, i, i, ratio);
    }
    
    // Convert V to GSL format
    gsl_matrix* V_gsl = gsl_matrix_calloc(p, p);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            gsl_matrix_set(V_gsl, i, j, V(i, j));
        }
    }
    
    // Calculate V * (Dsq/var) * V^T
    gsl_matrix* VDsqVarInv = gsl_matrix_calloc(p, p);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_gsl, DsqVarInv, 0.0, VDsqVarInv);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, VDsqVarInv, V_gsl, 0.0, XtOmegaX_full);
    
    gsl_matrix_free(DsqVarInv);
    gsl_matrix_free(VDsqVarInv);
    gsl_matrix_free(V_gsl);
  }
  
  NumericVector log10_BF(m);
  NumericVector log10_prior(m);
  NumericVector log10_posterior_score(m);
  NumericMatrix reg_weights;

  if (twas_weight) {
    reg_weights = NumericMatrix(m, p);
  }


  if (ss == 0 || ss == 1) {
    
    MLR mlr;
    mlr.init(yty, GtG, Gty, n, phi2_mat);

    for (int i = 0; i < m; i++) {
      // Extract the row
      vector<int> indicator = cmfg_matrix[i];
  
      vector<double> rst = mlr.compute_log10_BF(indicator, twas_weight);
      log10_BF[i] = rst[0];
  
      if (twas_weight) {
        for (int j = 0; j < p; j++) {
          reg_weights(i, j) = rst[j + 1];
        }
      }
  
      log10_prior[i] = compute_log10_prior(indicator, pi_vec);
      log10_posterior_score[i] = log10_BF[i] + log10_prior[i];
    }

  } else if (ss == 2) {
    // *** CHANGED: Infinitesimal model mode - use pre-computed XtOmegaX with MLR_infinitesimal class ***
    MLR_infinitesimal mlr_inf;
    mlr_inf.init(XtOmegaX_full, XtOmegay_full, phi2_mat);  // *** NEW: Initialize with pre-computed matrices ***
    
    for (int i = 0; i < m; i++) {
      // Extract the row
      vector<int> indicator = cmfg_matrix[i];
  
      vector<double> rst = mlr_inf.compute_log10_BF(indicator, twas_weight);  // *** Uses optimized method ***
      log10_BF[i] = rst[0];
  
      if (twas_weight) {
        for (int j = 0; j < p; j++) {
          reg_weights(i, j) = rst[j + 1];
        }
      }
  
      log10_prior[i] = compute_log10_prior(indicator, pi_vec);
      log10_posterior_score[i] = log10_BF[i] + log10_prior[i];
    }
  }

  // Free allocated memory
  gsl_matrix_free(GtG);
  gsl_matrix_free(Gty);
  gsl_matrix_free(XtOmegaX_full);
  gsl_vector_free(XtOmegay_full);

  return List::create(Named("log10_BF") = log10_BF,
                      Named("log10_prior") = log10_prior,
                      Named("log10_posterior_score") = log10_posterior_score,
                      Named("reg_weights") = reg_weights);
}