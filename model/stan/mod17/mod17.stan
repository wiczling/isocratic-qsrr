functions {
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
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
}

transformed data {
  vector[nAnalytes] logP_centered = logPobs - 2.5;

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
  real<lower=0, upper=1> alpha;
  vector<lower=0.01>[2] sdpi;
  matrix[nAnalytes, 2] param;
}
transformed parameters {
  matrix[nAnalytes, 2] miu;
  vector[nObs] logkx;
  matrix[2, 2] Omega = quad_form_diag(rho, omega);
  matrix[2, 2] L_Omega= cholesky_decompose(Omega); 
  miu[,1] = logkwHat + beta[1] * logP_centered + fgrp * pilogkw;
  miu[,2] = S1Hat + beta[2] * logP_centered + fgrp * piS1;
  
  for (i in 1:nAnalytes) {
    logkx[start[i]:end[i]] = funlogki(param[i,1], param[i,2], S2Hat, fi[start[i]:end[i]]);
  }
}

model {

  matrix[nAnalytes,nAnalytes] L_K=cholesky_decompose(similarity_x * alpha + (1-alpha)*identity_matrix(nAnalytes));
 
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
  sigma ~ normal(0, 0.05);
  alpha ~ normal(0.5,0.25);
  
  target += matrix_normal_lpdf(to_vector(param) | to_vector(miu), L_K, L_Omega);

  if (run_estimation == 1) {
    logkobs ~ student_t(7, logkx, sigma);
  }

}
