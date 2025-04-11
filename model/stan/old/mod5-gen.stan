 // Stan code that is based on
 // Comparison of chromatographic stationary phases using a Bayesian-based multilevel model
 // Simplified version of stan code that uses pKa as a predictor
 // one sigma
 // pKas and alphas without etas
 // reparametrized
 

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
    pHmpKa = log(10) *(2.6 - (pKaw + alpha * fi));
    
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
  int<lower=1> nSteps;
  array[nSteps,2] int idxch;
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
  vector[nObs] logkobs;                 // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
  array[nAnalytes] vector[maxR] alphax;
  array[nAnalytes] vector[maxR] pKawx;

  vector[11] fisim = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]';
  
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
  real<lower=0> S2Hat;               // typical value of S2a (log10 scale)
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  array[nSteps] vector[2] param;

}

transformed parameters {
  cov_matrix[2] Omega;
  array[nAnalytes] vector[maxR + 1] miulogkw;
  array[nAnalytes] vector[maxR + 1] miuS1;
  array[nSteps] vector[2] miu;
  
  array[nAnalytes] vector[maxR + 1] logkwx;
  array[nAnalytes] vector[maxR + 1] S1x;
  
  vector[nObs] logkx;

  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  {
  for (i in 1 : nAnalytes) {
    miulogkw[i, :] = (logkwHat + beta[1] * logPobs[i])*[1,1,1,1]';
    miuS1[i, :] = (S1Hat + beta[2] * logPobs[i])*[1,1,1,1]';
  }
  
  for (d in 1 : nA0) {
     miulogkw[idxA0[d], 1:4] += dlogkwHat;
     miuS1[idxA0[d], 1:4] += dS1Hat;
     }
  for (d in 1 : nB0) {
     miulogkw[idxB0[d], 1:4] += dlogkwHat;
     miuS1[idxB0[d], 1:4] += dS1Hat;
     }
   
  for (d in 1 : nGroupsA) {
      miulogkw[idxGroupsA[d,1], (idxGroupsA[d,2]+1) : 4] += dlogkwHat;
      miuS1[idxGroupsA[d,1], (idxGroupsA[d,2]+1) : 4] += dS1Hat;
      }
  
  for (d in 1 : nGroupsB) {
      miulogkw[idxGroupsB[d,1], 1:idxGroupsB[d,2]] += dlogkwHat;
      miuS1[idxGroupsB[d,1], 1:idxGroupsB[d,2]] += dS1Hat;
      }
   
// convert to long format
   for (nc in 1 : nSteps) { 
   miu[nc,1]  =  miulogkw[idxch[nc,1],idxch[nc,2]];
   miu[nc,2]  =  miuS1[idxch[nc,1],idxch[nc,2]];
  }
  }
  
  for (i in 1 : nAnalytes) { 
   logkwx[i, : ]  = [0,0,0,0]';
   S1x[i, : ]     = [0,0,0,0]';
  }
  
// convert long format to matrix
  for (nc in 1 : nSteps) { 
   logkwx[idxch[nc,1],idxch[nc,2]]  =  param[nc,1];
      S1x[idxch[nc,1],idxch[nc,2]]  =  param[nc,2];
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
  
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
 
  for (nc in 1 : nSteps) {
  param[nc] ~ multi_normal(miu[nc,1:2], Omega);
  }
  
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
  logkobs ~ student_t(7,logkx, sigma);
  }
}


generated quantities {
  
  array[nSteps] vector[2] sparam_ind;
  array[nSteps] vector[2] seta_ind;
  array[nSteps] vector[2] sparam_pop;
  array[nAnalytes] vector[maxR + 1] logkwx_pop;
  array[nAnalytes] vector[maxR + 1] S1x_pop;
  
  vector[nObs] slogkHat_ind;
  vector[nObs] slogkHat_pop;
  vector[nObs] slogk_ind;
  vector[nObs] slogk_pop;
  real<lower=0> sS2Hat;
  
  matrix[nAnalytes,11] simlogkHat_ind;
  matrix[nAnalytes,11] simlogkHat_pop;
  matrix[nAnalytes,11] simlogk_ind;
  matrix[nAnalytes,11] simlogk_pop;
  
  slogkHat_ind = logkx;
  sS2Hat = S2Hat;
  
  for (nc in 1 : nSteps) {
    sparam_pop[nc] = multi_normal_rng(miu[nc,1:2],Omega);
    sparam_ind[nc] = param[nc];
    seta_ind[nc] = (sparam_ind[nc] - miu[nc,1:2])./omega[1:2];
  }

  for (i in 1 : nAnalytes) { 
   logkwx_pop[i, : ]  = [0,0,0,0]';
   S1x_pop[i, : ]     = [0,0,0,0]';
  }
  
  for (nc in 1 : nSteps) { 
   logkwx_pop[idxch[nc,1],idxch[nc,2]]  =  sparam_pop[nc,1];
      S1x_pop[idxch[nc,1],idxch[nc,2]]  =  sparam_pop[nc,2];
  }
  
  for (z in 1 : nObs) {
   slogkHat_pop[z] = funlogki(logkwx_pop[analyte[z],:], S1x_pop[analyte[z],:], S2Hat, 
                      pKawx[analyte[z],:], alphax[analyte[z],:], R[analyte[z]],
                      fi[z]);
    slogk_ind[z] = student_t_rng(7, slogkHat_ind[z], sigma);
    slogk_pop[z] = student_t_rng(7, slogkHat_pop[z], sigma);
   }
   
  // dense grid simulations 
   for (f in 1 : 11) {
       for (a in 1 : nAnalytes) {
     simlogkHat_ind[a,f] = funlogki(logkwx[a,:], S1x[a,:], S2Hat, 
                      pKawx[a,:], alphax[a,:], R[a],
                      fisim[f]);
     simlogkHat_pop[a,f] = funlogki(logkwx_pop[a,:], S1x_pop[a,:], S2Hat, 
                      pKawx[a,:], alphax[a,:], R[a],
                      fisim[f]);
                      
     simlogk_ind[a,f] = student_t_rng(7, simlogkHat_ind[a,f], sigma);
     simlogk_pop[a,f] = student_t_rng(7, simlogkHat_pop[a,f], sigma);
       
       }}
   
}