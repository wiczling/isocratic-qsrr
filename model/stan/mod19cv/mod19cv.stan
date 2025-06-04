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


real matrix_normal_prec_lpdf(vector x, vector mu, matrix L_Prec_K, matrix L_Omega) {
    int n = rows(L_Prec_K); 
    int m = rows(L_Omega);  

    vector[n * m] z = x - mu;
    matrix[n, m] Z = to_matrix(z, n, m);

    // Apply L_Prec_K to rows instead of solving
    matrix[n, m] W = L_Prec_K * Z;              
    matrix[m, n] V = mdivide_left_tri_low(L_Omega, W'); // Solve L_Omega * V = W'
    real quad_form_x = dot_self(to_vector(V));

    real log_det_Prec_K = 2 * sum(log(diagonal(L_Prec_K)));
    real log_det_Omega = 2 * sum(log(diagonal(L_Omega)));

    // Since you're using precision matrix, flip sign of determinant contribution
    real log_det = -m * log_det_Prec_K + n * log_det_Omega;

    real constant = -0.5 * n * m * log(2.0 * pi());

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
   int nObscv; // subset of analytes and measurments
  array[nObscv] int idxcv; // indices
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
}

transformed data {
  vector[nAnalytes] logP_centered = logPobs - 2.5;
    vector[nObscv] logkobscv = logkobs[idxcv];
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

array[nAnalytes_uncorr] vector[2] param_uncorr;

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

  param[idx_corr[idx_corr_9],1:2]=param_corr_9;
  miu_corr_9 = miu[idx_corr[idx_corr_9], 1:2];

  param[idx_corr[idx_corr_10],1:2]=param_corr_10;
  miu_corr_10 = miu[idx_corr[idx_corr_10], 1:2];
  
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
 target += matrix_normal_prec_lpdf(to_vector(param_corr_1) | to_vector(miu_corr_1), L_K_prec_1, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_2 | point_mu_lower_2, point_sd_lower_2);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_2) | to_vector(miu_corr_2), L_K_prec_2, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_3 | point_mu_lower_3, point_sd_lower_3);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_3) | to_vector(miu_corr_3), L_K_prec_3, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_4 | point_mu_lower_4, point_sd_lower_4);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_4) | to_vector(miu_corr_4), L_K_prec_4, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_5 | point_mu_lower_5, point_sd_lower_5);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_5) | to_vector(miu_corr_5), L_K_prec_5, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_6 | point_mu_lower_6, point_sd_lower_6);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_6) | to_vector(miu_corr_6), L_K_prec_6, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_7 | point_mu_lower_7, point_sd_lower_7);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_7) | to_vector(miu_corr_7), L_K_prec_7, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_8 | point_mu_lower_8, point_sd_lower_8);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_8) | to_vector(miu_corr_8), L_K_prec_8, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_9 | point_mu_lower_9, point_sd_lower_9);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_9) | to_vector(miu_corr_9), L_K_prec_9, L_Omega);

  target += lkj_corr_cholesky_prec_point_lower_tri_lpdf(L_K_prec_10 | point_mu_lower_10, point_sd_lower_10);
 target += matrix_normal_prec_lpdf(to_vector(param_corr_10) | to_vector(miu_corr_10), L_K_prec_10, L_Omega);

 target += multi_normal_cholesky_lpdf(param_uncorr | miu_uncorr, L_Omega);

  if (run_estimation == 1) {
    target += student_t_lpdf(logkobscv | 7, logkx[idxcv],sigma);
  }

  
  // priors
  target += lprior;
}
