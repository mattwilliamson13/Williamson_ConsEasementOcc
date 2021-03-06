functions{
/**
* Return log probability of a unit-scale proper conditional autoregressive 
  * (CAR) prior with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int nobs, int W_n) {
      row_vector[nobs] phit_D; // phi' * D
      row_vector[nobs] phit_W; // phi' * W
      vector[nobs] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, nobs);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] += phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] += phi[W_sparse[i, 1]];
      }
    
      for (i in 1:nobs) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (sum(ldet_terms)
                    - (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data {
  // site-level occupancy covariates
  int<lower = 1> n_site;
  int<lower = 1> m_psi;
  matrix[n_site, m_psi] X_tct;
  
  // survey-level detection covariates
  int<lower = 1> total_surveys;
  int<lower = 1> m_p;
  matrix[total_surveys, m_p] X_bg;

  // survey level information  
  int<lower = 1, upper = n_site> site[total_surveys];
  int<lower = 0, upper = 1> y[total_surveys];
  int<lower = 0, upper = total_surveys> start_idx[n_site];
  int<lower = 0, upper = total_surveys> end_idx[n_site];
  
  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[n_site];
  
  // number of surveys at each site
  int<lower = 0> n_survey[n_site];
  
  //adjacency data
  matrix<lower = 0, upper = 1>[total_surveys, total_surveys] W_bg; //adjacency matrix bg
  int W_n_bg; //number of adjacent pairs bg
  //real<lower = 0, upper =1> alpha_occ;
  //real<lower = 0, upper =1> alpha_det;

}
transformed data {
  int W_sparse_det[W_n_bg, 2];
  vector[total_surveys] D_sparse_det;
  vector[total_surveys] lambda_det;
  { // generate sparse representation for W_det
  int counter_det;
  counter_det = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(total_surveys - 1)) {
      for (j in (i + 1):total_surveys) {
        if (W_bg[i, j] == 1) {
          W_sparse_det[counter_det, 1] = i;
          W_sparse_det[counter_det, 2] = j;
          counter_det += 1;
        }
      }
    }
  }
  for (i in 1:total_surveys) D_sparse_det[i] = sum(W_bg[i]);
  {
    vector[total_surveys] invsqrtD_det;  
    for (i in 1:total_surveys) {
      invsqrtD_det[i] = 1 / sqrt(D_sparse_det[i]);
    }
    lambda_det = eigenvalues_sym(quad_form(W_bg, diag_matrix(invsqrtD_det)));
  }
}

parameters {
  vector[total_surveys] phi_det;
  real<lower = 0, upper =0.999> alpha_det;
  real<lower = 0> sigma_det;

  vector[m_psi] beta_psi;
  vector[m_p] beta_p;
}
transformed parameters {
  vector[total_surveys] logit_p = X_bg * beta_p + phi_det * sigma_det;
  vector[n_site] logit_psi = X_tct * beta_psi;
}
model {
  vector[n_site] log_psi = log_inv_logit(logit_psi);
  vector[n_site] log1m_psi = log1m_inv_logit(logit_psi);

  phi_det ~ sparse_car(alpha_det, W_sparse_det, D_sparse_det, lambda_det, total_surveys, W_n_bg);
  alpha_det ~ beta(4,1);
  sigma_det ~ normal(0, 1.5);
  
  beta_psi ~ student_t(7.763,0, 1.566);
  beta_p ~ student_t(7.763,0, 1.566);
  
  for (i in 1:n_site) {
    if (n_survey[i] > 0) {
      if (any_seen[i]) {
        // site is occupied
        target += log_psi[i] 
                  + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] | 
                                         logit_p[start_idx[i]:end_idx[i]]);
      } else {
        // site may or may not be occupied
        target += log_sum_exp(
          log_psi[i] + bernoulli_logit_lpmf(y[start_idx[i]:end_idx[i]] |
                                            logit_p[start_idx[i]:end_idx[i]]), 
          log1m_psi[i]
        );
      }
    }
  }
}
