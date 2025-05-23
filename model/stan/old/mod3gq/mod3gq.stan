// Stan code to model the isocratic data
// Simplified version
 
functions {
// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower,
                                     real point_scale_lower) {
// works only for [2x2 matrix]
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
  int nOrphan; 
  array[nObs] int analyte;
  array[nAnalytes] int parents;
  vector[nObs] fi; 
  array[nAnalytes] int<lower = 1> start;
  array[nAnalytes] int<lower = 1> end;
  vector[nAnalytes] logPobs;
  int<lower=0> nK;             // number of predictors (functional groups)
  matrix[nAnalytes, nK] fgrp;     // predictor matrix (functional groups)
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
}

parameters {
  real logkwHat;            // typical logkw
  real S1Hat;               // effect of ACN on logkw
  real<lower=0> S2Hat;      // typical value of S2
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega1; // sd of BAV [logkw,S1]
  vector<lower=0>[2] omega2; // sd of BAV [logkw,S1]
  corr_matrix[2] rho1;       // correlation matrix [logkw vs. S1] 
  corr_matrix[2] rho2;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  vector[nK] pilogkw;        // pi-logkw
  vector[nK] piS1;           // pi-S1
  vector<lower = 0.01>[2] sdpi;     // between analyte variabilities for fgrp
  array[nAnalytes] vector[2] param;
}

transformed parameters {
  cov_matrix[2] Omega1;
  cov_matrix[2] Omega2;
  array[nAnalytes] vector[2] miu;
  vector[nObs] logkx;
  
  Omega1 = quad_form_diag(rho1, omega1); // diag_matrix(omega) * rho * diag_matrix(omega)
  Omega2 = quad_form_diag(rho2, omega2); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1 : nOrphan) {
  miu[i, 1] = logkwHat + beta[1] * (logPobs[i]-2.5)+ fgrp[i,1:nK] * pilogkw;
  miu[i, 2] = S1Hat    + beta[2] * (logPobs[i]-2.5)+ fgrp[i,1:nK] * piS1;
  }
  
  for (i in (nOrphan+1) : nAnalytes) {
  miu[i, 1] = logkwHat + beta[1] * (logPobs[i]-2.5)+ fgrp[i,1:nK] * pilogkw - miu[parents[i],1] ;
  miu[i, 2] = S1Hat    + beta[2] * (logPobs[i]-2.5)+ fgrp[i,1:nK] * piS1 - miu[parents[i],2] ;
  }
  
  for (i in 1 : nAnalytes) {
  logkx[start[i]:end[i]] = funlogki(param[i,1],param[i,2], S2Hat, fi[start[i]:end[i]]);
  }
}

model {
      
  logkwHat ~ normal(2, 4);
  S1Hat ~ normal(4, 2);
  S2Hat ~ lognormal(0.693, 0.125);
  beta[{1}] ~ normal(0.7, 0.125);
  beta[{2}] ~ normal(0.5, 0.5);
  omega1 ~ normal(0, 1);
  omega2 ~ normal(0, 1);
  rho1 ~ lkj_corr_point_lower_tri(0.75, 0.125);
  rho2 ~ lkj_corr_point_lower_tri(0, 0.25);
  
  pilogkw ~ normal(0,sdpi[1]);
  piS1    ~ normal(0,sdpi[2]);
  sdpi    ~ normal(0,0.5);
  
  for (i in 1 : nOrphan) {
  param[i] ~ multi_normal(miu[i], Omega1);
  }
  
  for (i in (nOrphan+1) : nAnalytes) {
  param[i] ~ multi_normal(param[parents[i]]+ miu[i], Omega2);
  }
  
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
   logkobs ~ student_t(7,logkx, sigma);
  }
}

generated quantities {
  
  array[nAnalytes] vector[2] sparam_ind;
  array[nAnalytes] vector[2] seta_ind;
  array[nAnalytes] vector[2] sparam_pop;
  vector[nObs] slogkHat_ind;
  vector[nObs] slogkHat_pop;
  vector[nObs] slogk_ind;
  vector[nObs] slogk_pop;
  real<lower=0> sS2Hat;
  
  sS2Hat = S2Hat;
  
   for (i in 1 : nOrphan) {
    sparam_pop[i] = multi_normal_rng(miu[i],Omega1);
    sparam_ind[i] = param[i];
    seta_ind[i] = param[i] - miu[i];
   }
  
   for (i in (nOrphan+1) : nAnalytes) {
   sparam_pop[i] = multi_normal_rng((sparam_pop[parents[i]] + miu[i]), Omega2);
   sparam_ind[i] = param[i];
   seta_ind[i] = param[i] - param[parents[i]] - miu[i];
  }
  
  for (i in 1 : nAnalytes) {
    slogkHat_ind[start[i]:end[i]] = funlogki(sparam_ind[i,1],sparam_ind[i,2], S2Hat, fi[start[i]:end[i]]);
    slogkHat_pop[start[i]:end[i]] = funlogki(sparam_pop[i,1],sparam_pop[i,2], S2Hat, fi[start[i]:end[i]]);
  }
  
   for (z in 1 : nObs) {
    slogk_ind[z] = student_t_rng(7, slogkHat_ind[z], sigma);
    slogk_pop[z] = student_t_rng(7, slogkHat_pop[z], sigma);
   }
}

