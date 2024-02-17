// Stan code to model the isocratic data
// Include logP and pKa as predictors. 
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
  array[nObs] int analyte;
  vector[nObs] fi; 
  array[nAnalytes] int<lower = 1> start;
  array[nAnalytes] int<lower = 1> end;
  vector[nAnalytes] logPobs;
  array[nAnalytes] vector[3] pKaslit;
  array[nAnalytes] vector[4] charges;
  int<lower=0> K;                       // number of predictors (functional groups)
  matrix[nAnalytes, K] fgrp;     // predictor matrix (functional groups) 
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
    array[nAnalytes] real fr;
    vector[4] cpHmpKa;
    vector[3] pHmpKa;
    array[nAnalytes] vector[4] log_lambda;
    
    for (i in 1 : nAnalytes) {
     pHmpKa = log(10) *(2.66-pKaslit[i]);
     cpHmpKa[1]=0;
     cpHmpKa[2:4]=cumulative_sum(pHmpKa);
     log_lambda[i] = log_softmax(cpHmpKa);
    }
}

parameters {
  real logkwHat;            // typical logkw
  real S1Hat;               // effect of ACN on logkw
  real dlogkwHat;           // typical dlogkw
  real dS1Hat;              // effect of ACN on dlogkw
  real<lower=0> S2Hat;      // typical value of S2
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  vector<lower=0>[2] kappa; // sd of BAV [dlogkw,dS1]
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  vector[K] pilogkw;        // pi-logkw
  vector[K] piS1;           // pi-S1
  vector<lower = 0.01>[2] sdpi;     // between analyte variabilities for fgrp
  array[nAnalytes] vector[2] param;
  vector[nAnalytes] etadlogkw;
  vector[nAnalytes] etadS1;
}

transformed parameters {
  cov_matrix[2] Omega;
  array[nAnalytes] vector[2] miu;
  vector[nObs] logkx;
  real logki;
  real S1i;
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
   for (i in 1 : nAnalytes) {
     
  logki = log_sum_exp(log_lambda[i] + log(10)*(logkwHat + (dlogkwHat + kappa[1]*etadlogkw[i]).*charges[i]))/log(10);
  S1i = log_sum_exp(log_lambda[i] + log(10)*(S1Hat + (dS1Hat + kappa[2]*etadS1[i]).*charges[i]))/log(10);
  
    miu[i, 1] = logki + beta[1] * logPobs[i] + fgrp[i,1:K] * pilogkw;
    miu[i, 2] = S1i   + beta[2] * logPobs[i] + fgrp[i,1:K] * piS1;
  }
  
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
  kappa ~ normal(0, 0.25);
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
 
  for (i in 1 : nAnalytes) {
  param[i] ~ multi_normal(miu[i,1:2], Omega);
  }
  
  etadlogkw ~ normal(0,1);
  etadS1 ~ normal(0,1);
  
  pilogkw ~ normal(0,sdpi[1]);
  piS1    ~ normal(0,sdpi[2]);
  sdpi ~ normal(0,0.1);
    
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
  array[nAnalytes] vector[2] miu_pop;
  real logki_pop;
  real S1i_pop;
  vector[nAnalytes] sdlogkw_pop;
  vector[nAnalytes] sdS1_pop;
 
  sS2Hat = S2Hat;
  
  for (i in 1 : nAnalytes) {
     
  sdlogkw_pop[i] = normal_rng(dlogkwHat,kappa[1]);
  sdS1_pop[i] = normal_rng(dS1Hat,kappa[2]);
  
  logki_pop = log_sum_exp(log_lambda[i] + log(10)*(logkwHat + sdlogkw_pop[i].*charges[i]))/log(10);
  S1i_pop = log_sum_exp(log_lambda[i] + log(10)*(S1Hat + sdS1_pop[i].*charges[i]))/log(10);
  
  miu_pop[i, 1] = logki_pop + beta[1] * logPobs[i] + fgrp[i,1:K] * pilogkw;
  miu_pop[i, 2] = S1i_pop   + beta[2] * logPobs[i] + fgrp[i,1:K] * piS1;
  }
  
  for(i in 1:nAnalytes){
    sparam_pop[i] = multi_normal_rng(miu_pop[i,1:2],Omega);
    sparam_ind[i] = param[i];
    seta_ind[i] = (sparam_ind[i] - miu[i,1:2])./omega[1:2];
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

