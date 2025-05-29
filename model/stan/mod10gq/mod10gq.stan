functions {

 vector lower_tri(matrix mat) {
    int d = rows(mat);
    int lower_tri_d = d * (d - 1) / 2;
    vector[lower_tri_d] lowerx;
    int countx = 1;

    for (c in 1:(d - 1)) {
      for (r in (c + 1):d) {
        lowerx[countx] = mat[r, c];
        countx += 1;
      }
    }

    return lowerx;
  }

real lkj_corr_cholesky_prec_point_lower_tri_lpdf(matrix cor_L, vector point_mu_lower, vector point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(lower_tri(cor) | point_mu_lower, point_scale_lower);
    return(lpdf);
 }

real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }

real kronecker_normal_lpdf(vector x, vector mu, matrix L_K_prec, matrix L_Omega) {
  int n = rows(L_K_prec);
  int m = rows(L_Omega);
  vector[n * m] z = x - mu;
  matrix[n, m] Z = to_matrix(z, n, m);
  matrix[n, m] W = mdivide_left_tri_low(L_K_prec, Z); // L_K_prec * W = Z
  matrix[m, n] V = mdivide_left_tri_low(L_Omega, W'); // L_Omega * V = W'
  real quad_form_x = sum(square(to_vector(V)));
  real log_det = -2 * m * sum(log(diagonal(L_K_prec))) + n * sum(log(diagonal(L_Omega)));
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
  int<lower=1> nD;
  vector[nD] point_mu_lower;
  vector[nD] point_sd_lower;
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
  vector<lower=0.01>[2] sdpi;
  matrix[nAnalytes_corr, 2] param_corr;
  array[nAnalytes_uncorr] vector[2] param_uncorr;
  real<lower=2> nu;
  cholesky_factor_cov[2] L_Omega_W;
  cholesky_factor_corr[nAnalytes_corr] L_K_prec;
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
  lprior += lognormal_lpdf(S2Hat | 0.693, 0.125);
  lprior += normal_lpdf(beta | [0.7, 0.5]', [0.125, 0.5]');
  lprior += normal_lpdf(omega | 0, 1);
  lprior += normal_lpdf(pilogkw | 0, sdpi[1]);
  lprior += normal_lpdf(piS1 | 0, sdpi[2]);
  lprior += normal_lpdf(sdpi | 0, 0.1);
  lprior += lkj_corr_point_lower_tri_lpdf(rho | 0.75, 0.125);
  lprior += gamma_lpdf(nu | 2, 0.1);
  lprior += normal_lpdf(sigma| 0, 0.05);

  real nprior=0;
  nprior += inv_wishart_cholesky_lpdf(L_Omega_W | nu, L_Omega);
  
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
  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec | point_mu_lower, point_sd_lower);
  target += kronecker_normal_lpdf(to_vector(param_corr) | to_vector(miu_corr), L_K_prec, L_Omega_W);
  target += multi_normal_cholesky_lpdf(param_uncorr | miu_uncorr, L_Omega_W);
  
  if (run_estimation == 1) {
    target += student_t_lpdf(logkobs | 7, logkx, sigma);
  }
  
  // latent parameters 
  target += nprior;
  
  // priors
  target += lprior;
}
generated quantities {
  
  matrix[nAnalytes,2] seta_ind;
  matrix[nAnalytes,2] sparam_pop;
  matrix[nAnalytes,2] sparam_ind;
  vector[nObs] slogkHat_ind;
  vector[nObs] slogkHat_pop;
  vector[nObs] slogk_ind;
  vector[nObs] slogk_pop;
  real<lower=0> sS2Hat;
  matrix[nAnalytes_uncorr,2] param_uncorr_pop;
  matrix[nAnalytes_corr, 2] param_corr_pop;  
  sS2Hat = S2Hat;
  sparam_ind = param;

 {
  matrix[nAnalytes_corr, 2] eta_pop;
  matrix[2, 2] rhoc;
   
for (i in 1 : nAnalytes_corr){
    for (j in 1 : 2) {
    eta_pop[i,j]=normal_rng(0,1);}}

rhoc= cholesky_decompose(rho);

matrix[nAnalytes_corr,nAnalytes_corr] L_K = mdivide_right_tri_low(identity_matrix(nAnalytes_corr), L_K_prec)';
param_corr_pop= miu_corr + L_K * eta_pop * diag_pre_multiply(omega, rhoc)';
}


 for (i in 1 : nAnalytes_uncorr) {
  param_uncorr_pop[i,1:2] = multi_normal_rng(miu[idx_uncorr[i],1:2], Omega)';
  }
  
  sparam_pop[idx_corr,1:2]=param_corr_pop;
  sparam_pop[idx_uncorr,1:2]=param_uncorr_pop;
  
  seta_ind = sparam_ind - miu;

  for (i in 1 : nAnalytes) {
    slogkHat_ind[start[i]:end[i]] = funlogki(sparam_ind[i,1],sparam_ind[i,2], S2Hat, fi[start[i]:end[i]]);
    slogkHat_pop[start[i]:end[i]] = funlogki(sparam_pop[i,1],sparam_pop[i,2], S2Hat, fi[start[i]:end[i]]);
  }
  
   for (z in 1 : nObs) {
    slogk_ind[z] = student_t_rng(7, slogkHat_ind[z], sigma);
    slogk_pop[z] = student_t_rng(7, slogkHat_pop[z], sigma);
   }
}