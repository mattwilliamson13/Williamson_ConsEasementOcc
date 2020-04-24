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
  
  //cty-level identifiers for varying intercept
  int<lower = 1> n_cty;
  int <lower = 1, upper = n_cty> tctincty[n_site];

  // summary of whether species is known to be present at each site
  int<lower = 0, upper = 1> any_seen[n_site];
  
  // number of surveys at each site
  int<lower = 0> n_survey[n_site];
  
  matrix<lower = 0, upper = 1>[n_site, n_site] W_tct; //adjacency matrix tract
  int W_n_tct; //number of adjacent pairs
}
transformed data {
  int W_sparse_occ[W_n_tct, 2];   // adjacency pairs
  vector[n_site] D_sparse_occ;     // diagonal of D (number of neigbors for each site)

  vector[n_site] lambda_occ;       // eigenvalues of invsqrtD * W * invsqrtD
  { // generate sparse representation for W
  int counter_occ;
  counter_occ = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n_site - 1)) {
      for (j in (i + 1):n_site) {
        if (W_tct[i, j] == 1) {
          W_sparse_occ[counter_occ, 1] = i;
          W_sparse_occ[counter_occ, 2] = j;
          counter_occ += 1;
        }
      }
    }
  }
  for (i in 1:n_site) D_sparse_occ[i] = sum(W_tct[i]);
  {
    vector[n_site] invsqrtD_occ;  
    for (i in 1:n_site) {
      invsqrtD_occ[i] = 1 / sqrt(D_sparse_occ[i]);
    }
    lambda_occ = eigenvalues_sym(quad_form(W_tct, diag_matrix(invsqrtD_occ)));
  }
  
}
parameters {
  vector[n_site] phi_occ;
  real<lower = 0, upper =0.999> alpha_occ;
  real<lower = 0> sigma_occ;
  vector[m_psi] beta_psi;
  real mu_alpha;
  vector[n_cty] alpha_cty_tilde;
  real<lower = 0> sigma_cty;
}
transformed parameters {
  vector[n_cty] alpha_cty = mu_alpha + sigma_cty * alpha_cty_tilde;
  vector[n_site] logit_psi = alpha_cty[tctincty] + X_tct * beta_psi + phi_occ * sigma_occ;
}
model {
  phi_occ ~ sparse_car(alpha_occ, W_sparse_occ, D_sparse_occ, lambda_occ, n_site, W_n_tct);
  alpha_occ ~ beta(3,2);
  sigma_occ ~ normal(0, 1.5);
  beta_psi ~ student_t(7.763,0, 1.566);
  mu_alpha ~ normal(0, 1.25);
  alpha_cty_tilde ~ std_normal();
  sigma_cty ~ student_t(4,0,1);
  any_seen ~ binomial_logit(n_survey, logit_psi);
}
