 // Stan code that is based on
 // Comparison of chromatographic stationary phases using a Bayesian-based multilevel model
 // Simplified version of stan code that uses pKa as a predictor
 // one sigma
 // pKas and alphas without etas

functions {

 // credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
  real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower,
                                     real point_scale_lower) {
    // works only for [2x2 matrix]
    real lpdf = lkj_corr_lpdf(rho | 1)
                + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }
  
 real funlogki(vector logkw, vector S1, real S2, vector pKaw, vector alpha, int R, real fi) {
    real lnki;
    real logki;
    vector[4] logkix;
    vector[3] pHmpKa;
    
    logkix = log(10) *(logkw - S1*(1+S2) * fi / (1 + S2 * fi));
    pHmpKa = log(10) *(2.66 - (pKaw + alpha * fi));
    
    if (R == 0) {
      lnki = logkix[1];
    } else if (R == 1) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1])-
             log1p_exp(pHmpKa[1]);
    } else if (R == 2) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1] + 
             log1p_exp(pHmpKa[2]+logkix[3]-logkix[2]))-
             log1p_exp(pHmpKa[1] + 
             log1p_exp(pHmpKa[2]));
    } else if (R == 3) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1] + 
             log1p_exp(pHmpKa[2]+logkix[3]-logkix[2] +
             log1p_exp(pHmpKa[3]+logkix[4]-logkix[3])))-
             log1p_exp(pHmpKa[1] + 
             log1p_exp(pHmpKa[2] +
             log1p_exp(pHmpKa[3])));
    }
    
    logki = lnki/log(10);
    
    return logki;
  } 

}


data {
  int nAnalytes;        
  int nObs;                 
  array[nObs] int analyte;
  array[nObs] real fi; 
  vector[nAnalytes] logPobs;
  int<lower=0, upper=3> maxR;
  array[nAnalytes] int<lower=0, upper=3> R;
  int<lower=1> nGroupsA;
  int<lower=1> nGroupsB;
  vector[nGroupsA] pKaslitA;
  vector[nGroupsB] pKaslitB;
  array[nGroupsA,2] int idxGroupsA;
  array[nGroupsB,2] int idxGroupsB;
  int<lower=1> nA0;
  int<lower=1> nB0;
  array[nA0] int idxA0;
  array[nB0] int idxB0; 
  int<lower=0> K;                       // number of predictors (functional groups)
  matrix[nAnalytes, K] fgrp;     // predictor matrix (functional groups) 
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
  array[nAnalytes] vector[maxR] alphax;
  array[nAnalytes] vector[maxR] pKawx;
  int grainsize = 1;
  
  for (i in 1 : nAnalytes) { 
   pKawx[i, : ]  = [0,0,0]';
   alphax[i, : ] = [0,0,0]'; 
  }
  
  for (d in 1 : nGroupsA) {
  pKawx[idxGroupsA[d,1], idxGroupsA[d,2]] = pKaslitA[d];
  alphax[idxGroupsA[d,1], idxGroupsA[d,2]]= 2;
  }
  
  for (d in 1 : nGroupsB) {
  pKawx[idxGroupsB[d,1], idxGroupsB[d,2]] = pKaslitB[d];
  alphax[idxGroupsB[d,1], idxGroupsB[d,2]]= -1;
}
  
  array[nObs] int ind = rep_array(1, nObs);
  
}
parameters {
  real logkwHat;            // typical logkw [Neutral]
  real S1Hat;               // effect of ACN on logkw [Neutral]
  real dlogkwHat;  // effect of dissociation on logkw [Acids, Bases]
  real dS1Hat;     // effect of dissociation on S1 [Acids, Bases] 
  real S2Hat;            // typical value of S2a (log10 scale)
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  vector<lower=0>[2] kappa; // sd of BAV [dlogkw,dS1]
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  vector[K] pilogkw;        // pi-logkw
  vector[K] piS1;           // pi-S1
  vector<lower = 0.01>[2] sdpi;     // between analyte variabilities for fgrp
  array[nAnalytes] vector[2] paramN;
  vector[nA0] etadlogkwA0;
  vector[nA0] etadS1A0;
  vector[nB0] etadlogkwB0;
  vector[nB0] etadS1B0;
  vector[nGroupsA] etadlogkwA;
  vector[nGroupsA] etadS1A;
  vector[nGroupsB] etadlogkwB;
  vector[nGroupsB] etadS1B;
}

transformed parameters {
  cov_matrix[2] Omega;
  array[nAnalytes] vector[2] miu;
  array[nAnalytes] vector[maxR + 1] logkwx;
  array[nAnalytes] vector[maxR + 1] S1x;
  vector[nObs] logkx;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + beta[1] * logPobs[i] + fgrp[i,1:K] * pilogkw;
    miu[i, 2] = S1Hat    + beta[2] * logPobs[i] + fgrp[i,1:K] * piS1;
  }
  for (i in 1 : nAnalytes) { 
   logkwx[i, : ]  =  paramN[i,1]*[1,1,1,1]';
      S1x[i, : ]  =  paramN[i,2]*[1,1,1,1]';
  }

  for (d in 1 : nA0) {
     logkwx[idxA0[d], 1:4] += dlogkwHat + kappa[1]*etadlogkwA0[d];
     S1x[idxA0[d], 1:4] += dS1Hat + kappa[2]*etadS1A0[d];
   }
    for (d in 1 : nB0) {
     logkwx[idxB0[d], 1:4] += dlogkwHat + kappa[1]*etadlogkwB0[d];
     S1x[idxB0[d], 1:4] += dS1Hat + kappa[2]*etadS1B0[d];
   }
   
  for (d in 1 : nGroupsA) {
      logkwx[idxGroupsA[d,1],  (idxGroupsA[d,2]+1) : 4] += dlogkwHat + kappa[1]*etadlogkwA[d];
         S1x[idxGroupsA[d,1],  (idxGroupsA[d,2]+1) : 4] += dS1Hat + kappa[2]*etadS1A[d];
         }
  
  for (d in 1 : nGroupsB) {
      logkwx[idxGroupsB[d,1], 1:idxGroupsB[d,2]] += dlogkwHat + kappa[1]*etadlogkwB[d];
      S1x[idxGroupsB[d,1],    1:idxGroupsB[d,2]] += dS1Hat + kappa[2]*etadS1B[d];
      }
  

  for(z in 1:nObs){
   logkx[z] = funlogki(logkwx[analyte[z],:],
   S1x[analyte[z],:], 
   S2Hat, 
   pKawx[analyte[z],:], 
   alphax[analyte[z],:], 
   R[analyte[z]],
   fi[z]);
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
 
  etadlogkwA0 ~ normal(0,1);
  etadS1A0 ~ normal(0,1);
  etadlogkwB0 ~ normal(0,1);
  etadS1B0 ~ normal(0,1);
  etadlogkwA ~ normal(0,1);
  etadS1A ~ normal(0,1);
  etadlogkwB ~ normal(0,1);
  etadS1B ~ normal(0,1);
  
  pilogkw ~ normal(0,sdpi[1]);
  piS1    ~ normal(0,sdpi[2]);
  sdpi ~ normal(0,0.1);
  
  for (i in 1 : nAnalytes) {
  paramN[i] ~ multi_normal(miu[i,1:2], Omega);
  }
  
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
  logkobs ~ student_t(7,logkx, sigma);
  }
}

generated quantities {
  
}
