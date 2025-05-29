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
  
  
real partial_sum(array[] int ind, int start, int end, 
vector logkobs_1,
vector fi_1,
array[] int starta_1,
array[] int enda_1,
vector point_mu_lower_1,
vector point_sd_lower_1,
matrix L_K_prec_1,
matrix param_corr_1,
matrix miu_corr_1,
 int nAnalyte_subset_1,
vector logkobs_2,
vector fi_2,
 array[] int starta_2,
array[] int enda_2,
vector point_mu_lower_2,
vector point_sd_lower_2,
matrix L_K_prec_2,
matrix param_corr_2,
matrix miu_corr_2,
 int nAnalyte_subset_2,
vector logkobs_3,
vector fi_3,
 array[] int starta_3,
array[] int enda_3,
vector point_mu_lower_3,
vector point_sd_lower_3,
matrix L_K_prec_3,
matrix param_corr_3,
matrix miu_corr_3,
 int nAnalyte_subset_3,
vector logkobs_4,
vector fi_4,
 array[] int starta_4,
array[] int enda_4,
vector point_mu_lower_4,
vector point_sd_lower_4,
matrix L_K_prec_4,
matrix param_corr_4,
matrix miu_corr_4,
 int nAnalyte_subset_4,
vector logkobs_5,
vector fi_5,
 array[] int starta_5,
array[] int enda_5,
vector point_mu_lower_5,
vector point_sd_lower_5,
matrix L_K_prec_5,
matrix param_corr_5,
matrix miu_corr_5,
 int nAnalyte_subset_5,
vector logkobs_6,
vector fi_6,
 array[] int starta_6,
array[] int enda_6,
vector point_mu_lower_6,
vector point_sd_lower_6,
matrix L_K_prec_6,
matrix param_corr_6,
matrix miu_corr_6,
 int nAnalyte_subset_6,
vector logkobs_7,
vector fi_7,
 array[] int starta_7,
array[] int enda_7,
vector point_mu_lower_7,
vector point_sd_lower_7,
matrix L_K_prec_7,
matrix param_corr_7,
matrix miu_corr_7,
 int nAnalyte_subset_7,
vector logkobs_8,
vector fi_8,
 array[] int starta_8,
array[] int enda_8,
vector point_mu_lower_8,
vector point_sd_lower_8,
matrix L_K_prec_8,
matrix param_corr_8,
matrix miu_corr_8,
 int nAnalyte_subset_8,
vector logkobs_9,
vector fi_9,
 array[] int starta_9,
array[] int enda_9,
vector point_mu_lower_9,
vector point_sd_lower_9,
matrix L_K_prec_9,
matrix param_corr_9,
matrix miu_corr_9,
int nAnalyte_subset_9,
vector logkobs_10,
vector fi_10,
array[] int starta_10,
array[] int enda_10,
vector point_mu_lower_10,
vector point_sd_lower_10,
matrix L_K_prec_10,
matrix param_corr_10,
matrix miu_corr_10,
 int nAnalyte_subset_10,
                   matrix L_Omega_W,
                   real S2Hat,
                   real sigma) {
                     
    real lp = 0;

for (z in start : end) {
  
if (z==1) {
int nobs = num_elements(logkobs_1);
vector[nobs] logkx_1;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_1| point_mu_lower_1, point_sd_lower_1);
lp += kronecker_normal_lpdf(to_vector(param_corr_1) | to_vector(miu_corr_1), L_K_prec_1, L_Omega_W);
for (i in 1:nAnalyte_subset_1) {
  logkx_1[starta_1[i]:enda_1[i]] =funlogki(param_corr_1[i,1], param_corr_1[i,2], S2Hat, fi_1[starta_1[i]:enda_1[i]]) ;
}
lp += student_t_lpdf(logkobs_1 | 7, logkx_1, sigma);
}
    
if (z==2) {
int nobs = num_elements(logkobs_2);
vector[nobs] logkx_2;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_2| point_mu_lower_2, point_sd_lower_2);
lp += kronecker_normal_lpdf(to_vector(param_corr_2) | to_vector(miu_corr_2), L_K_prec_2, L_Omega_W);
for (i in 1:nAnalyte_subset_2) {
  logkx_2[starta_2[i]:enda_2[i]] =funlogki(param_corr_2[i,1], param_corr_2[i,2], S2Hat, fi_2[starta_2[i]:enda_2[i]]) ;
}
lp += student_t_lpdf(logkobs_2 | 7, logkx_2, sigma);
}

if (z==3) {
int nobs = num_elements(logkobs_3);
vector[nobs] logkx_3;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_3| point_mu_lower_3, point_sd_lower_3);
lp += kronecker_normal_lpdf(to_vector(param_corr_3) | to_vector(miu_corr_3), L_K_prec_3, L_Omega_W);
for (i in 1:nAnalyte_subset_3) {
  logkx_3[starta_3[i]:enda_3[i]] =funlogki(param_corr_3[i,1], param_corr_3[i,2], S2Hat, fi_3[starta_3[i]:enda_3[i]]) ;
}
lp += student_t_lpdf(logkobs_3 | 7, logkx_3, sigma);
}

if (z==4) {
int nobs = num_elements(logkobs_4);
vector[nobs] logkx_4;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_4| point_mu_lower_4, point_sd_lower_4);
lp += kronecker_normal_lpdf(to_vector(param_corr_4) | to_vector(miu_corr_4), L_K_prec_4, L_Omega_W);
for (i in 1:nAnalyte_subset_4) {
  logkx_4[starta_4[i]:enda_4[i]] =funlogki(param_corr_4[i,1], param_corr_4[i,2], S2Hat, fi_4[starta_4[i]:enda_4[i]]) ;
}
lp += student_t_lpdf(logkobs_4 | 7, logkx_4, sigma);
}

if (z==5) {
  int nobs = num_elements(logkobs_5);
vector[nobs] logkx_5;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_5| point_mu_lower_5, point_sd_lower_5);
lp += kronecker_normal_lpdf(to_vector(param_corr_5) | to_vector(miu_corr_5), L_K_prec_5, L_Omega_W);
for (i in 1:nAnalyte_subset_5) {
  logkx_5[starta_5[i]:enda_5[i]] =funlogki(param_corr_5[i,1], param_corr_5[i,2], S2Hat, fi_5[starta_5[i]:enda_5[i]]) ;
}
lp += student_t_lpdf(logkobs_5 | 7, logkx_5, sigma);
}

if (z==6) {
  int nobs = num_elements(logkobs_6);
vector[nobs] logkx_6;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_6| point_mu_lower_6, point_sd_lower_6);
lp += kronecker_normal_lpdf(to_vector(param_corr_6) | to_vector(miu_corr_6), L_K_prec_6, L_Omega_W);
for (i in 1:nAnalyte_subset_6) {
  logkx_6[starta_6[i]:enda_6[i]] =funlogki(param_corr_6[i,1], param_corr_6[i,2], S2Hat, fi_6[starta_6[i]:enda_6[i]] );
}
lp += student_t_lpdf(logkobs_6 | 7, logkx_6, sigma);
}

if (z==7) {
  int nobs = num_elements(logkobs_7);
vector[nobs] logkx_7;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_7| point_mu_lower_7, point_sd_lower_7);
lp += kronecker_normal_lpdf(to_vector(param_corr_7) | to_vector(miu_corr_7), L_K_prec_7, L_Omega_W);
for (i in 1:nAnalyte_subset_7) {
  logkx_7[starta_7[i]:enda_7[i]] =funlogki(param_corr_7[i,1], param_corr_7[i,2], S2Hat, fi_7[starta_7[i]:enda_7[i]]) ;
}
lp += student_t_lpdf(logkobs_7 | 7, logkx_7, sigma);
}

if (z==8) {
  int nobs = num_elements(logkobs_8);
vector[nobs] logkx_8;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_8| point_mu_lower_8, point_sd_lower_8);
lp += kronecker_normal_lpdf(to_vector(param_corr_8) | to_vector(miu_corr_8), L_K_prec_8, L_Omega_W);
for (i in 1:nAnalyte_subset_8) {
  logkx_8[starta_8[i]:enda_8[i]] =funlogki(param_corr_8[i,1], param_corr_8[i,2], S2Hat, fi_8[starta_8[i]:enda_8[i]]) ;
}
lp += student_t_lpdf(logkobs_8 | 7, logkx_8, sigma);
}

if (z==9) {
  int nobs = num_elements(logkobs_9);
vector[nobs] logkx_9;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_9| point_mu_lower_9, point_sd_lower_9);
lp += kronecker_normal_lpdf(to_vector(param_corr_9) | to_vector(miu_corr_9), L_K_prec_9, L_Omega_W);
for (i in 1:nAnalyte_subset_9) {
  logkx_9[starta_9[i]:enda_9[i]] =funlogki(param_corr_9[i,1], param_corr_9[i,2], S2Hat, fi_9[starta_9[i]:enda_9[i]]) ;
}
lp += student_t_lpdf(logkobs_9 | 7, logkx_9, sigma);
}

if (z==10) {
  int nobs = num_elements(logkobs_10);
vector[nobs] logkx_10;
lp += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_10| point_mu_lower_10, point_sd_lower_10);
lp += kronecker_normal_lpdf(to_vector(param_corr_10) | to_vector(miu_corr_10), L_K_prec_10, L_Omega_W);
for (i in 1:nAnalyte_subset_10) {
  logkx_10[starta_10[i]:enda_10[i]] =funlogki(param_corr_10[i,1], param_corr_10[i,2], S2Hat, fi_10[starta_10[i]:enda_10[i]]);
}
lp += student_t_lpdf(logkobs_10 | 7, logkx_10, sigma);
}

return lp;
  }
}
}
data {
  int nAnalytes;
  int nObs;
  vector[nAnalytes] logPobs;
  int<lower=0> nK;
  matrix[nAnalytes, nK] fgrp;
  
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

int<lower=1> nD_9;
 int<lower=1> nAnalyte_subset_9;
 vector[nD_9] point_mu_lower_9;
 vector[nD_9] point_sd_lower_9;
 array[nAnalyte_subset_9] int idx_corr_9;

int<lower=1> nD_10;
 int<lower=1> nAnalyte_subset_10;
 vector[nD_10] point_mu_lower_10;
 vector[nD_10] point_sd_lower_10;
 array[nAnalyte_subset_10] int idx_corr_10;
 
int<lower=1> nAnalyte_subset_11; // uncorelated
array[nAnalyte_subset_11] int idx_corr_11; 

int<lower=1> nObs_1;
vector[nObs_1] fi_1;
vector[nObs_1] logkobs_1;
array[nAnalyte_subset_1] int<lower=1> starta_1;
array[nAnalyte_subset_1] int<lower=1> enda_1;

int<lower=1> nObs_2;
vector[nObs_2] fi_2;
vector[nObs_2] logkobs_2;
array[nAnalyte_subset_2] int<lower=1> starta_2;
array[nAnalyte_subset_2] int<lower=1> enda_2;

int<lower=1> nObs_3;
vector[nObs_3] fi_3;
vector[nObs_3] logkobs_3;
array[nAnalyte_subset_3] int<lower=1> starta_3;
array[nAnalyte_subset_3] int<lower=1> enda_3;

int<lower=1> nObs_4;
vector[nObs_4] fi_4;
vector[nObs_4] logkobs_4;
array[nAnalyte_subset_4] int<lower=1> starta_4;
array[nAnalyte_subset_4] int<lower=1> enda_4;

int<lower=1> nObs_5;
vector[nObs_5] fi_5;
vector[nObs_5] logkobs_5;
array[nAnalyte_subset_5] int<lower=1> starta_5;
array[nAnalyte_subset_5] int<lower=1> enda_5;

int<lower=1> nObs_6;
vector[nObs_6] fi_6;
vector[nObs_6] logkobs_6;
array[nAnalyte_subset_6] int<lower=1> starta_6;
array[nAnalyte_subset_6] int<lower=1> enda_6;

int<lower=1> nObs_7;
vector[nObs_7] fi_7;
vector[nObs_7] logkobs_7;
array[nAnalyte_subset_7] int<lower=1> starta_7;
array[nAnalyte_subset_7] int<lower=1> enda_7;

int<lower=1> nObs_8;
vector[nObs_8] fi_8;
vector[nObs_8] logkobs_8;
array[nAnalyte_subset_8] int<lower=1> starta_8;
array[nAnalyte_subset_8] int<lower=1> enda_8;

int<lower=1> nObs_9;
vector[nObs_9] fi_9;
vector[nObs_9] logkobs_9;
array[nAnalyte_subset_9] int<lower=1> starta_9;
array[nAnalyte_subset_9] int<lower=1> enda_9;

int<lower=1> nObs_10;
vector[nObs_10] fi_10;
vector[nObs_10] logkobs_10;
array[nAnalyte_subset_10] int<lower=1> starta_10;
array[nAnalyte_subset_10] int<lower=1> enda_10;

int<lower=1> nObs_11;
vector[nObs_11] fi_11;
vector[nObs_11] logkobs_11;
array[nAnalyte_subset_11] int<lower=1> starta_11;
array[nAnalyte_subset_11] int<lower=1> enda_11;
}

transformed data {
  vector[nAnalytes] logP_centered = logPobs - 2.5;
  int grainsize = 1;
  array[10] int ind = rep_array(1, 10);
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
  matrix[nAnalyte_subset_9, 2] param_corr_9;
  matrix[nAnalyte_subset_10, 2] param_corr_10;
  array[nAnalyte_subset_11] vector[2] param_11;
  real<lower=2> nu;
  cholesky_factor_cov[2] L_Omega_W;
  cholesky_factor_corr[nAnalyte_subset_1] L_K_prec_1;
  cholesky_factor_corr[nAnalyte_subset_2] L_K_prec_2;
  cholesky_factor_corr[nAnalyte_subset_3] L_K_prec_3;
  cholesky_factor_corr[nAnalyte_subset_4] L_K_prec_4;
  cholesky_factor_corr[nAnalyte_subset_5] L_K_prec_5;
  cholesky_factor_corr[nAnalyte_subset_6] L_K_prec_6;
  cholesky_factor_corr[nAnalyte_subset_7] L_K_prec_7;
  cholesky_factor_corr[nAnalyte_subset_8] L_K_prec_8;
  cholesky_factor_corr[nAnalyte_subset_9] L_K_prec_9;
  cholesky_factor_corr[nAnalyte_subset_10] L_K_prec_10;

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
matrix[nAnalyte_subset_9, 2] miu_corr_9;
matrix[nAnalyte_subset_10, 2] miu_corr_10;
array[nAnalyte_subset_11] vector[2] miu_11;
vector[nObs_11] logkx_11;

  matrix[2, 2] Omega = quad_form_diag(rho, omega);
  cholesky_factor_cov[2] L_Omega= cholesky_decompose(Omega); 
  
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
  
 miu_corr_1 =  miu[idx_corr_1, 1:2];
 miu_corr_2 = miu[idx_corr_2, 1:2];
 miu_corr_3 = miu[idx_corr_3, 1:2];
 miu_corr_4 = miu[idx_corr_4, 1:2];
 miu_corr_5 = miu[idx_corr_5, 1:2];
 miu_corr_6 = miu[idx_corr_6, 1:2];
 miu_corr_7 = miu[idx_corr_7, 1:2];
 miu_corr_8 = miu[idx_corr_8, 1:2];
 miu_corr_9 = miu[idx_corr_9, 1:2];
 miu_corr_10 = miu[idx_corr_10, 1:2];
 for (i in 1:nAnalyte_subset_11) {miu_11[i] = miu[idx_corr_11[i], 1:2]';}
  
  for (i in 1:nAnalyte_subset_11) {
    logkx_11[starta_11[i]:enda_11[i]] = funlogki(param_11[i,1], param_11[i,2], S2Hat, fi_11[starta_11[i]:enda_11[i]]);
  }

}

model {
  
  //target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_1 | point_mu_lower_1, point_sd_lower_1);
  //target += kronecker_normal_lpdf(to_vector(param_corr_1) | to_vector(miu_corr_1), L_K_prec_1, L_Omega_W);
  
target += reduce_sum(partial_sum, ind, grainsize, 
 logkobs_1,
 fi_1,
 starta_1,
 enda_1,
 point_mu_lower_1,
 point_sd_lower_1,
 L_K_prec_1,
 param_corr_1,
 miu_corr_1,
 nAnalyte_subset_1,
 logkobs_2,
 fi_2,
 starta_2,
 enda_2,
 point_mu_lower_2,
 point_sd_lower_2,
 L_K_prec_2,
 param_corr_2,
 miu_corr_2,
 nAnalyte_subset_2,
 logkobs_3,
 fi_3,
 starta_3,
 enda_3,
 point_mu_lower_3,
 point_sd_lower_3,
 L_K_prec_3,
 param_corr_3,
 miu_corr_3,
 nAnalyte_subset_3,
 logkobs_4,
 fi_4,
 starta_4,
 enda_4,
 point_mu_lower_4,
 point_sd_lower_4,
 L_K_prec_4,
 param_corr_4,
 miu_corr_4,
 nAnalyte_subset_4,
 logkobs_5,
 fi_5,
 starta_5,
 enda_5,
 point_mu_lower_5,
 point_sd_lower_5,
 L_K_prec_5,
 param_corr_5,
 miu_corr_5,
 nAnalyte_subset_5,
 logkobs_6,
 fi_6,
 starta_6,
 enda_6,
 point_mu_lower_6,
 point_sd_lower_6,
 L_K_prec_6,
 param_corr_6,
 miu_corr_6,
 nAnalyte_subset_6,
 logkobs_7,
 fi_7,
 starta_7,
 enda_7,
 point_mu_lower_7,
 point_sd_lower_7,
 L_K_prec_7,
 param_corr_7,
 miu_corr_7,
 nAnalyte_subset_7,
 logkobs_8,
 fi_8,
 starta_8,
 enda_8,
 point_mu_lower_8,
 point_sd_lower_8,
 L_K_prec_8,
 param_corr_8,
 miu_corr_8,
 nAnalyte_subset_8,
 logkobs_9,
 fi_9,
 starta_9,
 enda_9,
 point_mu_lower_9,
 point_sd_lower_9,
 L_K_prec_9,
 param_corr_9,
 miu_corr_9,
 nAnalyte_subset_9,
 logkobs_10,
 fi_10,
 starta_10,
 enda_10,
 point_mu_lower_10,
 point_sd_lower_10,
 L_K_prec_10,
 param_corr_10,
 miu_corr_10,
 nAnalyte_subset_10,
 L_Omega_W,
 S2Hat,
 sigma);
                         
  target += multi_normal_cholesky_lpdf(param_11 | miu_11, L_Omega_W);
  target += student_t_lpdf(logkobs_11 | 7, logkx_11, sigma);

  // latent parameters 
  target += nprior;
  
  // priors
  target += lprior;
}
