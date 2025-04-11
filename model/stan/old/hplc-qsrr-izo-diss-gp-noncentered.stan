// Stan code to model the isocratic data
// Include logP and pKa as predictors. 
// Simplified version
 
functions {
// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
  real lkj_corr_cholesky_point_lower_tri_lpdf(matrix rhoc, real point_mu_lower, real point_scale_lower) {
// works only for [2x2 matrix]

    matrix[2,2] rho = multiply_lower_tri_self_transpose(rhoc);

    real lpdf = lkj_corr_lpdf(rho | 1)
                + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }
  

 vector funlogki(real logkw, real S1, real S2, vector fi) {
// isocratic Neue model
    int d = num_elements(fi);
    vector[d] logki;
    logki = logkw-S1*(1+S2)*fi./(1+S2*fi);
    return logki;
  }
}

data {
  int nAnalytes;        
  int nObs;                 
  array[nObs] int analyte;
  vector[nObs] fi; 
  array[nAnalytes] int<lower = 1> start;
  array[nAnalytes] int<lower = 1> end;
  vector[nAnalytes] logPobs;
  array[nAnalytes] vector[3] pKaslit;
  array[nAnalytes] vector[4] charges;
  vector[nObs] logkobs;                 // observed retention factors 
  matrix[nAnalytes, nAnalytes] distance_square; // distance 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
    array[nAnalytes] real fr;
    vector[4] cpHmpKa;
    vector[4] lambda;
    vector[3] pHmpKa;
    
    real delta =0.001;  
    for (i in 1 : nAnalytes) {
     pHmpKa = log(10) *(2.66-pKaslit[i]);
     cpHmpKa[1]=0;
     cpHmpKa[2:4]=cumulative_sum(pHmpKa);
     lambda = softmax(cpHmpKa);
     fr[i] = sum(lambda.*charges[i]);
    }

}

parameters {
  real logkwHat;            // typical logkw
  real S1Hat;               // effect of ACN on logkw
  real dlogkwHat;            // typical dlogkw
  real dS1Hat;               // effect of ACN on dlogkw
  real<lower=0> S2Hat;      // typical value of S2
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  cholesky_factor_corr[2] rhoc;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  matrix[nAnalytes,2] eta;
    real<lower=0, upper=1> rho_gp;
}

transformed parameters {

  matrix[nAnalytes,2] param;
  matrix[nAnalytes,2] miu;
  vector[nObs] logkx;
  matrix[nAnalytes, nAnalytes] K;
  matrix[nAnalytes, nAnalytes] L_K;
  
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + fr[i]*dlogkwHat + beta[1] * logPobs[i];
    miu[i, 2] = S1Hat    + fr[i]*dS1Hat + beta[2] * logPobs[i];
  }
  
  {
   for (i in 1:(nAnalytes - 1)) {
     K[i, i] = 1 + delta;
   for (j in (i + 1):nAnalytes) {
      K[i, j] = exp(-0.5/square(rho_gp)*distance_square[i,j]);
      K[j, i] = K[i, j];
     }
   }
  
  K[nAnalytes, nAnalytes] = 1 + delta;
  L_K = cholesky_decompose(K);
  }
  
  param = miu + (L_K * eta * diag_pre_multiply(omega, rhoc)');
     
  for (i in 1 : nAnalytes) {
    logkx[start[i]:end[i]] = funlogki(param[i,1],param[i,2], S2Hat, fi[start[i]:end[i]]);
  }

}

model {
  logkwHat ~ normal(2, 4);
  S1Hat ~ normal(4, 2);
  dlogkwHat ~ normal(-1, 0.25);
  dS1Hat ~ normal(0, 0.5);
  S2Hat ~ lognormal(0.693, 0.125);
  beta[{1}] ~ normal(0.7, 0.125);
  beta[{2}] ~ normal(0.5, 0.5);
  omega ~ normal(0, 1);
  
  rhoc ~ lkj_corr_cholesky_point_lower_tri(0.75, 0.125);
 
  //to_vector(param) ~ multi_normal(to_vector(miu), OmegaK);
  //to_vector(param) ~  multi_normal_cholesky(to_vector(miu), L_K);
  to_vector(eta) ~ std_normal();
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
   logkobs ~ student_t(7,logkx, sigma);
  }
}

generated quantities {
  
  corr_matrix[2] rho;
  rho = rhoc * rhoc';
}
