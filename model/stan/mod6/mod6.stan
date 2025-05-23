functions {
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }

real kronecker_normal_lpdf(vector x, vector mu, matrix L_K, matrix L_Omega) {
  int n = rows(L_K);
  int m = rows(L_Omega);
  vector[n * m] z = x - mu;
  matrix[n, m] Z = to_matrix(z, n, m);
  matrix[n, m] W = mdivide_left_tri_low(L_K, Z); // L_K * W = Z
  matrix[m, n] V = mdivide_left_tri_low(L_Omega, W'); // L_Omega * V = W'
  real quad_form_x = sum(square(to_vector(V)));
  real log_det = m * sum(log(diagonal(L_K))) + n * sum(log(diagonal(L_Omega)));
  return -0.5 * (quad_form_x + log_det);
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
  matrix[nAnalytes,nAnalytes] distance_x;
}

transformed data {
  vector[nAnalytes] logP_centered = logPobs - 2.5;
  array[nAnalytes_corr] vector[nAnalytes_corr] distance_s;
  
  for (i in 1:nAnalytes_corr) {
  distance_s[i] = distance_x[idx_corr,idx_corr[i]];
 }
}

parameters {
  real logkwHat;
  real S1Hat;
  real<lower=0> S2Hat;
  vector[2] beta;
  vector<lower=0>[2] omega;
  corr_matrix[2] rho;
  real<lower=0> sigma;
  vector[nK] pilogkw;
  vector[nK] piS1;
  vector<lower=0.01>[2] sdpi;
  matrix[nAnalytes_corr, 2] param_corr;
  matrix[nAnalytes_uncorr, 2] param_uncorr;
  real<lower=0> rhogp;
}

transformed parameters {
  matrix[nAnalytes, 2] miu;
  vector[nObs] logkx;
  matrix[2, 2] Omega = quad_form_diag(rho, omega);
  matrix[2, 2] L_Omega= cholesky_decompose(Omega); 
  matrix[nAnalytes, 2] param;
  miu[,1] = logkwHat + beta[1] * logP_centered + fgrp * pilogkw;
  miu[,2] = S1Hat + beta[2] * logP_centered + fgrp * piS1;
  
  param[idx_corr,1:2]=param_corr;
  param[idx_uncorr,1:2]=param_uncorr;
  
  for (i in 1:nAnalytes) {
    logkx[start[i]:end[i]] = funlogki(param[i,1], param[i,2], S2Hat, fi[start[i]:end[i]]);
  }
}

model {
  matrix[nAnalytes_corr, nAnalytes_corr] L_K;
  
  {
  matrix[nAnalytes_corr, nAnalytes_corr] K = gp_exp_quad_cov(distance_s, 1.0, rhogp);
  for (n in 1:nAnalytes_corr) {
    K[n, n] = K[n, n] + 0.01;
  }
  L_K = cholesky_decompose(K);
  }
 
  
  logkwHat ~ normal(2, 4);
  S1Hat ~ normal(4, 2);
  S2Hat ~ lognormal(0.693, 0.125);
  beta[1] ~ normal(0.7, 0.125);
  beta[2] ~ normal(0.5, 0.5);
  omega ~ normal(0, 1);
  pilogkw ~ normal(0, sdpi[1]);
  piS1 ~ normal(0, sdpi[2]);
  sdpi ~ normal(0, 0.1);
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
  
  target += kronecker_normal_lpdf(to_vector(param_corr) | to_vector(miu[idx_corr,1:2]), L_K, L_Omega);

  for (i in 1 : nAnalytes_uncorr) {
  param_uncorr[i] ~ multi_normal(miu[idx_uncorr[i],1:2], Omega);
  }
  
  sigma ~ normal(0, 0.05);
  if (run_estimation == 1) {
    logkobs ~ student_t(7, logkx, sigma);
  }
  
  rhogp ~ inv_gamma(5, 2);
}
