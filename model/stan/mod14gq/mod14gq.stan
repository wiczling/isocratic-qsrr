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
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 2);
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
  
 int<lower=1> nD_1;
 int<lower=1> nAnalyte_subset_1;
 vector[nD_1] point_mu_lower_1;
 vector[nD_1] point_sd_lower_1;
 array[nAnalyte_subset_1] int idx_corr_1;

int<lower=1> nD_2;
 int<lower=1> nAnalyte_subset_2;
 vector[nD_2] point_mu_lower_2;
 vector[nD_2] point_sd_lower_2;
 array[nAnalyte_subset_2] int idx_corr_2;

int<lower=1> nD_3;
 int<lower=1> nAnalyte_subset_3;
 vector[nD_3] point_mu_lower_3;
 vector[nD_3] point_sd_lower_3;
 array[nAnalyte_subset_3] int idx_corr_3;

int<lower=1> nD_4;
 int<lower=1> nAnalyte_subset_4;
 vector[nD_4] point_mu_lower_4;
 vector[nD_4] point_sd_lower_4;
 array[nAnalyte_subset_4] int idx_corr_4;

int<lower=1> nD_5;
 int<lower=1> nAnalyte_subset_5;
 vector[nD_5] point_mu_lower_5;
 vector[nD_5] point_sd_lower_5;
 array[nAnalyte_subset_5] int idx_corr_5;

int<lower=1> nD_6;
 int<lower=1> nAnalyte_subset_6;
 vector[nD_6] point_mu_lower_6;
 vector[nD_6] point_sd_lower_6;
 array[nAnalyte_subset_6] int idx_corr_6;

int<lower=1> nD_7;
 int<lower=1> nAnalyte_subset_7;
 vector[nD_7] point_mu_lower_7;
 vector[nD_7] point_sd_lower_7;
 array[nAnalyte_subset_7] int idx_corr_7;

int<lower=1> nD_8;
 int<lower=1> nAnalyte_subset_8;
 vector[nD_8] point_mu_lower_8;
 vector[nD_8] point_sd_lower_8;
 array[nAnalyte_subset_8] int idx_corr_8;

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
  matrix[nAnalyte_subset_1, 2] param_corr_1;
matrix[nAnalyte_subset_2, 2] param_corr_2;
matrix[nAnalyte_subset_3, 2] param_corr_3;
matrix[nAnalyte_subset_4, 2] param_corr_4;
matrix[nAnalyte_subset_5, 2] param_corr_5;
matrix[nAnalyte_subset_6, 2] param_corr_6;
matrix[nAnalyte_subset_7, 2] param_corr_7;
matrix[nAnalyte_subset_8, 2] param_corr_8;

array[nAnalytes_uncorr] vector[2] param_uncorr;
cholesky_factor_corr[nAnalyte_subset_1] L_K_prec_1;
cholesky_factor_corr[nAnalyte_subset_2] L_K_prec_2;
cholesky_factor_corr[nAnalyte_subset_3] L_K_prec_3;
cholesky_factor_corr[nAnalyte_subset_4] L_K_prec_4;
cholesky_factor_corr[nAnalyte_subset_5] L_K_prec_5;
cholesky_factor_corr[nAnalyte_subset_6] L_K_prec_6;
cholesky_factor_corr[nAnalyte_subset_7] L_K_prec_7;
cholesky_factor_corr[nAnalyte_subset_8] L_K_prec_8;

}

transformed parameters {
  matrix[nAnalytes, 2] miu;
matrix[nAnalyte_subset_1, 2] miu_corr_1;
matrix[nAnalyte_subset_2, 2] miu_corr_2;
matrix[nAnalyte_subset_3, 2] miu_corr_3;
matrix[nAnalyte_subset_4, 2] miu_corr_4;
matrix[nAnalyte_subset_5, 2] miu_corr_5;
matrix[nAnalyte_subset_6, 2] miu_corr_6;
matrix[nAnalyte_subset_7, 2] miu_corr_7;
matrix[nAnalyte_subset_8, 2] miu_corr_8;

  
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
  lprior += normal_lpdf(sigma| 0, 0.05);
  
  miu[,1] = logkwHat + beta[1] * logP_centered + fgrp * pilogkw;
  miu[,2] = S1Hat + beta[2] * logP_centered + fgrp * piS1;
  
  param[idx_corr[idx_corr_1],1:2]=param_corr_1;
  miu_corr_1 = miu[idx_corr[idx_corr_1], 1:2];

  param[idx_corr[idx_corr_2],1:2]=param_corr_2;
  miu_corr_2 = miu[idx_corr[idx_corr_2], 1:2];

  param[idx_corr[idx_corr_3],1:2]=param_corr_3;
  miu_corr_3 = miu[idx_corr[idx_corr_3], 1:2];

  param[idx_corr[idx_corr_4],1:2]=param_corr_4;
  miu_corr_4 = miu[idx_corr[idx_corr_4], 1:2];

  param[idx_corr[idx_corr_5],1:2]=param_corr_5;
  miu_corr_5 = miu[idx_corr[idx_corr_5], 1:2];

  param[idx_corr[idx_corr_6],1:2]=param_corr_6;
  miu_corr_6 = miu[idx_corr[idx_corr_6], 1:2];

  param[idx_corr[idx_corr_7],1:2]=param_corr_7;
  miu_corr_7 = miu[idx_corr[idx_corr_7], 1:2];

  param[idx_corr[idx_corr_8],1:2]=param_corr_8;
  miu_corr_8 = miu[idx_corr[idx_corr_8], 1:2];
  
  for (i in 1:nAnalytes_uncorr) {
  param[idx_uncorr[i], 1:2] = param_uncorr[i]';
  miu_uncorr[i] = miu[idx_uncorr[i], 1:2]';
  }
  
  for (i in 1:nAnalytes) {
    logkx[start[i]:end[i]] = funlogki(param[i,1], param[i,2], S2Hat, fi[start[i]:end[i]]);
  }

}

model {
  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_1 | point_mu_lower_1, point_sd_lower_1);
 target += kronecker_normal_lpdf(to_vector(param_corr_1) | to_vector(miu_corr_1), L_K_prec_1, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_2 | point_mu_lower_2, point_sd_lower_2);
 target += kronecker_normal_lpdf(to_vector(param_corr_2) | to_vector(miu_corr_2), L_K_prec_2, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_3 | point_mu_lower_3, point_sd_lower_3);
 target += kronecker_normal_lpdf(to_vector(param_corr_3) | to_vector(miu_corr_3), L_K_prec_3, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_4 | point_mu_lower_4, point_sd_lower_4);
 target += kronecker_normal_lpdf(to_vector(param_corr_4) | to_vector(miu_corr_4), L_K_prec_4, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_5 | point_mu_lower_5, point_sd_lower_5);
 target += kronecker_normal_lpdf(to_vector(param_corr_5) | to_vector(miu_corr_5), L_K_prec_5, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_6 | point_mu_lower_6, point_sd_lower_6);
 target += kronecker_normal_lpdf(to_vector(param_corr_6) | to_vector(miu_corr_6), L_K_prec_6, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_7 | point_mu_lower_7, point_sd_lower_7);
 target += kronecker_normal_lpdf(to_vector(param_corr_7) | to_vector(miu_corr_7), L_K_prec_7, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_8 | point_mu_lower_8, point_sd_lower_8);
 target += kronecker_normal_lpdf(to_vector(param_corr_8) | to_vector(miu_corr_8), L_K_prec_8, L_Omega);

 target += multi_normal_cholesky_lpdf(param_uncorr | miu_uncorr, L_Omega);

  if (run_estimation == 1) {
    target += student_t_lpdf(logkobs | 7, logkx, sigma);
  }
  
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
  matrix[nAnalyte_subset_1, 2] param_corr_pop_1;
  matrix[nAnalyte_subset_2, 2] param_corr_pop_2;
 matrix[nAnalyte_subset_3, 2] param_corr_pop_3;
 matrix[nAnalyte_subset_4, 2] param_corr_pop_4;
 matrix[nAnalyte_subset_5, 2] param_corr_pop_5;
 matrix[nAnalyte_subset_6, 2] param_corr_pop_6;
 matrix[nAnalyte_subset_7, 2] param_corr_pop_7;
 matrix[nAnalyte_subset_8, 2] param_corr_pop_8;
  sS2Hat = S2Hat;
  sparam_ind = param;
matrix[nAnalyte_subset_1,nAnalyte_subset_1]  rho1 = L_K_prec_1 * L_K_prec_1';
matrix[nAnalyte_subset_2,nAnalyte_subset_2]  rho2 = L_K_prec_2 * L_K_prec_2';
matrix[nAnalyte_subset_3,nAnalyte_subset_3]  rho3 = L_K_prec_3 * L_K_prec_3';
matrix[nAnalyte_subset_4,nAnalyte_subset_4]  rho4 = L_K_prec_4 * L_K_prec_4';
matrix[nAnalyte_subset_5,nAnalyte_subset_5]  rho5 = L_K_prec_5 * L_K_prec_5';
matrix[nAnalyte_subset_6,nAnalyte_subset_6]  rho6 = L_K_prec_6 * L_K_prec_6';
matrix[nAnalyte_subset_7,nAnalyte_subset_7]  rho7 = L_K_prec_7 * L_K_prec_7';
matrix[nAnalyte_subset_8,nAnalyte_subset_8]  rho8 = L_K_prec_8 * L_K_prec_8';

 {
  matrix[nAnalyte_subset_1, 2] eta_pop_1;
  matrix[nAnalyte_subset_2, 2] eta_pop_2;
  matrix[nAnalyte_subset_3, 2] eta_pop_3;
  matrix[nAnalyte_subset_4, 2] eta_pop_4;
  matrix[nAnalyte_subset_5, 2] eta_pop_5;
  matrix[nAnalyte_subset_6, 2] eta_pop_6;
  matrix[nAnalyte_subset_7, 2] eta_pop_7;
  matrix[nAnalyte_subset_8, 2] eta_pop_8;
  matrix[2, 2] rhoc;
   
for (i in 1 : nAnalyte_subset_1){
    for (j in 1 : 2) {
    eta_pop_1[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_2){
    for (j in 1 : 2) {
    eta_pop_2[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_3){
    for (j in 1 : 2) {
    eta_pop_3[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_4){
    for (j in 1 : 2) {
    eta_pop_4[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_5){
    for (j in 1 : 2) {
    eta_pop_5[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_6){
    for (j in 1 : 2) {
    eta_pop_6[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_7){
    for (j in 1 : 2) {
    eta_pop_7[i,j]=normal_rng(0,1);}}
 for (i in 1 : nAnalyte_subset_8){
    for (j in 1 : 2) {
    eta_pop_8[i,j]=normal_rng(0,1);}}

rhoc= cholesky_decompose(rho);

matrix[nAnalyte_subset_1,nAnalyte_subset_1] L_K_1 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_1), L_K_prec_1)';
matrix[nAnalyte_subset_2,nAnalyte_subset_2] L_K_2 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_2), L_K_prec_2)';
matrix[nAnalyte_subset_3,nAnalyte_subset_3] L_K_3 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_3), L_K_prec_3)';
matrix[nAnalyte_subset_4,nAnalyte_subset_4] L_K_4 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_4), L_K_prec_4)';
matrix[nAnalyte_subset_5,nAnalyte_subset_5] L_K_5 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_5), L_K_prec_5)';
matrix[nAnalyte_subset_6,nAnalyte_subset_6] L_K_6 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_6), L_K_prec_6)';
matrix[nAnalyte_subset_7,nAnalyte_subset_7] L_K_7 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_7), L_K_prec_7)';
matrix[nAnalyte_subset_8,nAnalyte_subset_8] L_K_8 = mdivide_right_tri_low(identity_matrix(nAnalyte_subset_8), L_K_prec_8)';

param_corr_pop_1= miu_corr_1 + L_K_1  * eta_pop_1 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_2= miu_corr_2 + L_K_2  * eta_pop_2 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_3= miu_corr_3 + L_K_3  * eta_pop_3 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_4= miu_corr_4 + L_K_4  * eta_pop_4 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_5= miu_corr_5 + L_K_5  * eta_pop_5 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_6= miu_corr_6 + L_K_6  * eta_pop_6 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_7= miu_corr_7 + L_K_7  * eta_pop_7 * diag_pre_multiply(omega, rhoc)';
param_corr_pop_8= miu_corr_8 + L_K_8  * eta_pop_8 * diag_pre_multiply(omega, rhoc)';
}


 for (i in 1 : nAnalytes_uncorr) {
  param_uncorr_pop[i,1:2] = multi_normal_rng(miu[idx_uncorr[i],1:2], Omega)';
  }
  
sparam_pop[idx_corr[idx_corr_1],1:2]=param_corr_pop_1;
sparam_pop[idx_corr[idx_corr_2],1:2]=param_corr_pop_2;
sparam_pop[idx_corr[idx_corr_3],1:2]=param_corr_pop_3;
sparam_pop[idx_corr[idx_corr_4],1:2]=param_corr_pop_4;
sparam_pop[idx_corr[idx_corr_5],1:2]=param_corr_pop_5;
sparam_pop[idx_corr[idx_corr_6],1:2]=param_corr_pop_6;
sparam_pop[idx_corr[idx_corr_7],1:2]=param_corr_pop_7;
sparam_pop[idx_corr[idx_corr_8],1:2]=param_corr_pop_8;

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
