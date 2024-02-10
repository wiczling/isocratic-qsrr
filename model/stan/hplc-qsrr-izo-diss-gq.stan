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
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
    array[nAnalytes] real fr;
    vector[4] cpHmpKa;
    vector[4] lambda;
    vector[3] pHmpKa;
    
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
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  array[nAnalytes] vector[2] param;
}

transformed parameters {
  cov_matrix[2] Omega;
  array[nAnalytes] vector[2] miu;
  vector[nObs] logkx;
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + fr[i]*dlogkwHat + beta[1] * logPobs[i];
    miu[i, 2] = S1Hat    + fr[i]*dS1Hat + beta[2] * logPobs[i];
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
  
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
 
  for (i in 1 : nAnalytes) {
  param[i] ~ multi_normal(miu[i,1:2], Omega);
  }
  
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
   logkobs ~ student_t(7,logkx, sigma);
  }
}

generated quantities {
  
  array[nAnalytes] vector[2] param_pred;
  
  for(i in 1:nAnalytes){
    param_pred[i] = multi_normal_rng(miu[i,1:2],Omega); 
  }
  
}
