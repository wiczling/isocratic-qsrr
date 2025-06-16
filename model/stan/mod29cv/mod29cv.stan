functions {
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }

  real matrix_t_lpdf(matrix X,
                     matrix mu,
                     matrix L_SigmaR,
                     matrix L_SigmaC,
                     real nu) {
    
    real ldSigmaR = 2 * sum(log(diagonal(L_SigmaR)));
    real ldSigmaC = 2 * sum(log(diagonal(L_SigmaC)));
    
    int p = rows(X);
    int q = cols(X);
    int pq = p * q;
    
    real nuq = nu + q - 1.0;
    real nupq = nuq + p;
    
    matrix[p, q] Z = X - mu;
    matrix[p, q] A = mdivide_left_tri_low(L_SigmaR, Z); 
    
    matrix[q, q] SigmaC = multiply_lower_tri_self_transpose(L_SigmaC); 
    matrix[q, q] Bq_raw = SigmaC + Z' * A;
	matrix[q, q] Bq = 0.5 * (Bq_raw + Bq_raw');
    matrix[q, q] L_Bq = cholesky_decompose(Bq); 
    real ldBq = 2 * sum(log(diagonal(L_Bq)));
    
    real ldens = ldBq - ldSigmaC;
    ldens = 0.5 * nupq * ldens + 0.5 * q * ldSigmaR + 0.5 * p * ldSigmaC + 0.5 * pq * log(pi());
 
    return -ldens + lmgamma(q, 0.5 * nupq) - lmgamma(q, 0.5 * nuq);
  }
  
real matrix_normal_lpdf(vector x, vector mu, matrix L_K, matrix L_Omega) {
    int n = rows(L_K); 
    int m = rows(L_Omega);  

    vector[n * m] z = x - mu;
    matrix[n, m] Z = to_matrix(z, n, m);

    matrix[n, m] W = mdivide_left_tri_low(L_K, Z);     // Solve L_K * W = Z
    matrix[m, n] V = mdivide_left_tri_low(L_Omega, W'); // Solve L_Omega * V = W'
    real quad_form_x = dot_self(to_vector(V));


    real log_det_K = 2 * sum(log(diagonal(L_K)));
    real log_det_Omega = 2 * sum(log(diagonal(L_Omega)));
    real log_det = m * log_det_K + n * log_det_Omega;
    real constant = -0.5 * n * m * log(2.0 * pi());

    // Final log-likelihood
    return constant - 0.5 * (quad_form_x + log_det);
  }

  vector funlogki(real logkw, real S1, real S2, vector fi) {
    int d = num_elements(fi);
    vector[d] logki = logkw - S1 * (1 + S2) * fi ./ (1 + S2 * fi);
    return logki;
  }
}

data {
  int nAnalytes;
  int nAnalytes_corr;
  int nAnalytes_uncorr;
  array[nAnalytes_corr] int idx_corr;
  array[nAnalytes_uncorr] int idx_uncorr;
  int nfg_diss;
  int nfg_nondiss;
  array[nfg_diss] int idx_diss;
  array[nfg_nondiss] int idx_nondiss; 
  int nObs;
  array[nObs] int analyte;
  vector[nObs] fi;
  array[nAnalytes] int<lower=1> start;
  array[nAnalytes] int<lower=1> end;
  vector[nAnalytes] logPobs;
  int<lower=0> nK;
  matrix[nAnalytes, nK] fgrp;
  vector[nObs] logkobs;
  int<lower=0, upper=1> run_estimation;
  matrix[nAnalytes,nAnalytes] similarity_x;
  int mGroup;
  array[nAnalytes_corr] int<lower=1, upper=mGroup> group;
      int nObscv; // subset of analytes and measurments
  array[nObscv] int idxcv; // indices
}

transformed data {
  vector[nObscv] logkobscv = logkobs[idxcv];
  vector[nAnalytes] logP_centered = logPobs - 2.5;
  matrix[nAnalytes_corr,nAnalytes_corr] similarity_s = similarity_x[idx_corr,idx_corr]; 

  array[mGroup] int n_corr_per_group;
  array[mGroup, nAnalytes_corr] int idx_corr_group;

  for (g in 1:mGroup) {
  int count = 0;
  for (i in 1:nAnalytes_corr) {
    if (group[i] == g) {
      count += 1;
      idx_corr_group[g][count] = i;
    }
  }
  n_corr_per_group[g] = count;
}


}

parameters {
  real logkwHat;
  real S1Hat;
    real dS1Hat;
  real dlogkwHat;
  real<lower=0> S2Hat;
  vector[2] beta;
  vector<lower=0>[2] omega;
  corr_matrix[2] rho;
  real<lower=0> sigma;
  vector[nK] pilogkw;
  vector[nK] piS1;
  vector<lower=0.01>[4] sdpi;
  matrix[nAnalytes_corr, 2] param_corr;
  array[nAnalytes_uncorr] vector[2] param_uncorr;
  real<lower=0, upper=1> alpha;
  real<lower=2> nu;  
  cholesky_factor_cov[2] L_Omega_W;
}

transformed parameters {
  matrix[nAnalytes, 2] miu;
  matrix[nAnalytes_corr, 2] miu_corr;
  array[nAnalytes_uncorr] vector[2] miu_uncorr;
  vector[nObs] logkx;
  matrix[2, 2] Omega = quad_form_diag(rho, omega);
  cholesky_factor_cov[2] L_Omega= cholesky_decompose(Omega); 
  matrix[nAnalytes, 2] param;

  
  real lprior=0;
  lprior += normal_lpdf(logkwHat | 2, 4);
  lprior += normal_lpdf(S1Hat | 4, 2);
    lprior += normal_lpdf(dlogkwHat | -1, 0.5);
  lprior += normal_lpdf(dS1Hat | 0, 0.5);
  lprior += lognormal_lpdf(S2Hat | 0.693, 0.125);
  lprior += normal_lpdf(beta | [0.7, 0.5]', [0.125, 0.5]');
  lprior += normal_lpdf(omega | 0, 1);
  lprior += normal_lpdf(pilogkw[idx_nondiss] | 0, sdpi[1]);
  lprior += normal_lpdf(piS1[idx_nondiss] | 0, sdpi[2]);
  lprior += normal_lpdf(pilogkw[idx_diss] | dlogkwHat, sdpi[3]);
  lprior += normal_lpdf(piS1[idx_diss] | dS1Hat, sdpi[4]);
  lprior += normal_lpdf(sdpi | 0, 0.1);
  lprior += lkj_corr_point_lower_tri_lpdf(rho | 0.75, 0.125);
  lprior += normal_lpdf(alpha | 0.5, 0.25);
  lprior += normal_lpdf(sigma| 0, 0.05);

  
  miu[,1] = logkwHat + beta[1] * logP_centered + fgrp * pilogkw;
  miu[,2] = S1Hat + beta[2] * logP_centered + fgrp * piS1;
  
  param[idx_corr,1:2]=param_corr;
  for (i in 1:nAnalytes_uncorr) {
  param[idx_uncorr[i], 1:2] = param_uncorr[i]';
  }
  miu_corr = miu[idx_corr, 1:2];
   for (i in 1:nAnalytes_uncorr) {
  miu_uncorr[i] = miu[idx_uncorr[i], 1:2]';
  }
  
  for (i in 1:nAnalytes) {
    logkx[start[i]:end[i]] = funlogki(param[i,1], param[i,2], S2Hat, fi[start[i]:end[i]]);
  }

}

model {

  target += inv_wishart_cholesky_lpdf(L_Omega_W | nu, L_Omega);
  for (g in 1:mGroup) {
  int n_g = n_corr_per_group[g];
  array[n_g] int idx_g = idx_corr_group[g][1:n_g];
  matrix[n_g, 2] param_g;
  matrix[n_g, 2] miu_g;
  matrix[n_g, n_g] K_g;
  matrix[n_g, n_g] L_K_g;
  
  for (i in 1:n_g) {
    param_g[i] = param_corr[idx_g[i], 1:2];
    miu_g[i] = miu_corr[idx_g[i], 1:2];
  }

  L_K_g =cholesky_decompose(similarity_s[idx_g,idx_g]*alpha*nu + (1.0-alpha)*identity_matrix(n_g)*nu);
   target += matrix_normal_lpdf(to_vector(param_g) | to_vector(miu_g), L_K_g, L_Omega_W);
}

  target += multi_student_t_cholesky_lpdf(param_uncorr | nu, miu_uncorr, L_Omega);

  if (run_estimation == 1) {
    target += student_t_lpdf(logkobscv | 7, logkx[idxcv],sigma);
  }
  
  // prior 
  target += lprior;
}
