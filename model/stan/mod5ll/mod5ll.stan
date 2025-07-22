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
  array[nObs] int analyte;
  int nfg_diss;
  int nfg_nondiss;
  array[nfg_diss] int idx_diss;
  array[nfg_nondiss] int idx_nondiss; 
  vector[nObs] fi; 
  array[nAnalytes] int<lower = 1> start;
  array[nAnalytes] int<lower = 1> end;
  vector[nAnalytes] logPobs;
  int<lower=0> nK;                       // number of predictors (functional groups)
  matrix[nAnalytes, nK] fgrp;     // predictor matrix (functional groups) 
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  vector[nAnalytes] logP_centered = logPobs - 2.5;
}

parameters {
  real logkwHat;            // typical logkw
  real S1Hat;               // effect of ACN on logkw
  real dlogkwHat;
  real dS1Hat;
  real<lower=0> S2Hat;      // typical value of S2
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  vector[nK] pilogkw;        // pi-logkw
  vector[nK] piS1;           // pi-S1
  vector<lower = 0.01>[4] sdpi;     // between analyte variabilities for fgrp
  array[nAnalytes] vector[2] param;
}

transformed parameters {
  cov_matrix[2] Omega;
  matrix[nAnalytes, 2] miu;
  vector[nObs] logkx;
  
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
  lprior += normal_lpdf(sigma| 0, 0.05);
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  miu[,1] = logkwHat + beta[1] * logP_centered + fgrp * pilogkw;
  miu[,2] = S1Hat + beta[2] * logP_centered + fgrp * piS1;
  
  for (i in 1 : nAnalytes) {
  logkx[start[i]:end[i]] = funlogki(param[i,1],param[i,2], S2Hat, fi[start[i]:end[i]]);
  }

}

model {
  for (i in 1 : nAnalytes) {
   target += multi_normal_lpdf(param[i] | miu[i,1:2], Omega);
  }

  if (run_estimation == 1) {
    target += student_t_lpdf(logkobs | 7, logkx, sigma);
  }
  target += lprior;
}

generated quantities {
  vector[nObs] log_lik;
{
  vector[nObs] logkxx;
  for (i in 1 : nAnalytes) {
    logkxx[start[i]:end[i]] = funlogki(param[i,1],param[i,2], S2Hat, fi[start[i]:end[i]]);
  }
  
  for (z in 1 : nObs) {
   log_lik[z] = student_t_lpdf(logkobs[z] | 7, logkxx[z], sigma);  
  }
}
}

