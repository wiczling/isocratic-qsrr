# need to add charges for acids and bases that do not change
# chargesA[,1] and chargesB[,4]


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
             log1p_exp(pHmpKa[1] + logkix[2] - logkix[1]) -
             log1p_exp(pHmpKa[1]);
    } else if (R == 2) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1] + log1p_exp(pHmpKa[2]+logkix[3]-logkix[2])) -
             log1p_exp(pHmpKa[1] + log1p_exp(pHmpKa[2]));
    } else if (R == 3) {
      lnki = logkix[1] +
             log1p_exp(pHmpKa[1]+logkix[2]-logkix[1] + log1p_exp(pHmpKa[2]+logkix[3]-logkix[2] + log1p_exp(pHmpKa[3]+logkix[4]-logkix[3]))) -
             log1p_exp(pHmpKa[1] + log1p_exp(pHmpKa[2] + log1p_exp(pHmpKa[3])));
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
  int<lower=0, upper=3> maxR; //
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
  vector[nObs] logkobs;                   // observed retention factors 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}


parameters {
  
  real logkwHat; // typical logkw [Neutral]
  real S1Hat;   // effect of ACN on logkw [Neutral]
  array[2] real dlogkwHat; // effect of dissociation on logkw [Acids, Bases]
  array[2] real dS1Hat;   // effect of dissociation on S1 [Acids, Bases] 
  real logS2Hat; // typical value of S2a (log10 scale)
  vector[2] beta; // effect of logP on logkw and S1
  array[2] real alphaHat;   // effect of ACN on pKa [Acids, Bases]
  real<lower=0> tau;   // sd for between analyte variability of pKa's
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1m]
  corr_matrix[2] rho;      // correlation matrix [logkw vs. S1m] 
  vector<lower=0>[2] kappa; // sd of BAV [dlogkw,dS1a]
  
  array[nAnalytes] vector[2] paramN;
  vector[nGroupsA] dlogkwA;
  vector[nGroupsB] dlogkwB;
  vector[nGroupsA] dS1A;
  vector[nGroupsB] dS1B;
  
  vector[nA0] dlogkwA0;
  vector[nB0] dlogkwB0;
  vector[nA0] dS1A0;
  vector[nB0] dS1B0;
  
  // Dissociation
  vector[nGroupsA] pKawA;
  vector[nGroupsB] pKawB;

  // residual variability
  real<lower=0> msigma;   // typical sigma for the 1st column
  real<lower=0> ssigma;
  vector[nAnalytes] logsigma; 
}
transformed parameters {
  
  cov_matrix[2] Omega;
  array[nAnalytes] vector[2] miu;
  array[nAnalytes] vector[maxR + 1] logkwx;
  array[nAnalytes] vector[maxR + 1] S1x;
  real S2x;
  array[nAnalytes] vector[maxR] alphax;
  array[nAnalytes] vector[maxR] pKawx;
  vector[nObs] sigmax;
  vector[nObs] logkx;
  
  vector[nGroupsA] alphaA;
  vector[nGroupsB] alphaB;
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + beta[1] * (logPobs[i]-2.2);
    miu[i, 2] = S1Hat    + beta[2] * (logPobs[i]-2.2);
  }

  for (i in 1 : nAnalytes) { 
   logkwx[i, : ]  =  paramN[i,1]*[1,1,1,1]';
  }
  
   for (d in 1 : nA0) {
     logkwx[idxA0[d], 1] += dlogkwA0[d];
     logkwx[idxA0[d], 2] += dlogkwA0[d];
     logkwx[idxA0[d], 3] += dlogkwA0[d];
     logkwx[idxA0[d], 4] += dlogkwA0[d]; 
   }
    for (d in 1 : nB0) {
     logkwx[idxB0[d], 1] += dlogkwB0[d];
     logkwx[idxB0[d], 2] += dlogkwB0[d];
     logkwx[idxB0[d], 3] += dlogkwB0[d];
     logkwx[idxB0[d], 4] += dlogkwB0[d]; 
   }
   
  for (d in 1 : nGroupsA) {
    if (idxGroupsA[d,2]==1) {
      logkwx[idxGroupsA[d,1], 2] += dlogkwA[d];
      logkwx[idxGroupsA[d,1], 3] += dlogkwA[d];
      logkwx[idxGroupsA[d,1], 4] += dlogkwA[d];
   }
   if (idxGroupsA[d,2]==2) {
      logkwx[idxGroupsA[d,1], 3] += dlogkwA[d];
      logkwx[idxGroupsA[d,1], 4] += dlogkwA[d];
   }
   if (idxGroupsA[d,2]==3) {
      logkwx[idxGroupsA[d,1], 4] += dlogkwA[d];
   }}
  
  for (d in 1 : nGroupsB) {
      if (idxGroupsB[d,2]==3) {
      logkwx[idxGroupsB[d,1], 3] += dlogkwB[d];
      logkwx[idxGroupsB[d,1], 2] += dlogkwB[d];
      logkwx[idxGroupsB[d,1], 1] += dlogkwB[d];
   }
   if (idxGroupsB[d,2]==2) {
      logkwx[idxGroupsB[d,1], 2] += dlogkwB[d];
      logkwx[idxGroupsB[d,1], 1] += dlogkwB[d];
   }
   if (idxGroupsB[d,2]==1) {
      logkwx[idxGroupsB[d,1], 1] += dlogkwB[d];
   }}
    
    
  for (i in 1 : nAnalytes) { 
   S1x[i, : ]  =  paramN[i,2]*[1,1,1,1]';
   }
   
   for (d in 1 : nA0) {
     S1x[idxA0[d], 1] += dS1A0[d];
     S1x[idxA0[d], 2] += dS1A0[d];
     S1x[idxA0[d], 3] += dS1A0[d];
     S1x[idxA0[d], 4] += dS1A0[d]; 
   }
   
   for (d in 1 : nB0) {
     S1x[idxB0[d], 1] += dS1B0[d];
     S1x[idxB0[d], 2] += dS1B0[d];
     S1x[idxB0[d], 3] += dS1B0[d];
     S1x[idxB0[d], 4] += dS1B0[d]; 
   }
   
  for (d in 1 : nGroupsA) {
    if (idxGroupsA[d,2]==1) {
      S1x[idxGroupsA[d,1], 2] += dS1A[d];
      S1x[idxGroupsA[d,1], 3] += dS1A[d];
      S1x[idxGroupsA[d,1], 4] += dS1A[d];
   }
   if (S1x[d,2]==2) {
      S1x[idxGroupsA[d,1], 3] += dS1A[d];
      S1x[idxGroupsA[d,1], 4] += dS1A[d];
   }
   if (idxGroupsA[d,2]==3) {
      S1x[idxGroupsA[d,1], 4] += dS1A[d];
   }}
  
  for (d in 1 : nGroupsB) {
      if (idxGroupsB[d,2]==1) {
      S1x[idxGroupsB[d,1], 1] += dS1B[d];
      S1x[idxGroupsB[d,1], 2] += dS1B[d];
      S1x[idxGroupsB[d,1], 3] += dS1B[d];
   }
   if (idxGroupsB[d,2]==2) {
      S1x[idxGroupsB[d,1], 1] += dS1B[d];
      S1x[idxGroupsB[d,1], 2] += dS1B[d];
   }
   if (idxGroupsB[d,2]==3) {
      S1x[idxGroupsB[d,1], 1] += dS1B[d];
   }}
  
  S2x = 10^logS2Hat;

  for (i in 1 : nAnalytes) { 
   pKawx[i, : ]  = [0,0,0]';
   alphax[i, : ]  = [0,0,0]'; 
  }
  
  for (d in 1 : nGroupsA) {
  pKawx[idxGroupsA[d,1], idxGroupsA[d,2]] = pKawA[d];
  alphax[idxGroupsA[d,1], idxGroupsA[d,2]]= alphaHat[1];
  }
  
  for (d in 1 : nGroupsB) {
  pKawx[idxGroupsB[d,1], idxGroupsB[d,2]] = pKawB[d];
  alphax[idxGroupsB[d,1], idxGroupsB[d,2]]= alphaHat[2];
  }
   
  for (z in 1 : nObs) { 
   sigmax[z] = exp(logsigma[analyte[z]]);
   
   logkx[z] = funlogki(logkwx[analyte[z],:],
   S1x[analyte[z],:], 
   S2x, 
   pKawx[analyte[z],:], 
   alphax[analyte[z],:], 
   R[analyte[z]],
   fi[analyte[z]]);
  }
  
}

model {
  logkwHat ~ normal(2.2, 2);
  S1Hat ~ normal(4, 1);
  dlogkwHat ~ normal(-1, 0.125);
  dS1Hat   ~ normal(0, 0.5);
  logS2Hat ~ normal(0.3, 0.125);
  beta[{1}] ~ normal(1, 0.125);
  beta[{2}] ~ normal(0.5, 0.5);
  alphaHat[{1}] ~ normal(2, 0.25);
  alphaHat[{2}] ~ normal(-1, 0.25);
  tau ~ normal(0, 0.25);
  omega ~ normal(0, 2);
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
  kappa ~ normal(0, 0.25);
    
  for (i in 1 : nAnalytes) {
  paramN[i] ~ multi_normal(miu[i,1:2], Omega);
  }
  
  dlogkwA ~ normal(dlogkwHat[1], kappa[1]);
  dlogkwB ~ normal(dlogkwHat[2], kappa[1]);
  dS1A ~ normal(dS1Hat[1], kappa[2]);
  dS1B ~ normal(dS1Hat[2], kappa[2]);
  dlogkwA0 ~ normal(dlogkwHat[1], kappa[1]);
  dlogkwB0 ~ normal(dlogkwHat[2], kappa[1]);
  dS1A0 ~ normal(dS1Hat[1], kappa[2]);
  dS1B0 ~ normal(dS1Hat[2], kappa[2]);
  
  pKawA ~ normal(pKaslitA, tau);
  pKawB ~ normal(pKaslitB, tau);

  msigma ~ normal(0,0.1);
  ssigma ~ normal(0,0.2);
  logsigma  ~ normal(log(msigma),ssigma); 
  
  if (run_estimation == 1) {
   logkobs ~ student_t(7,logkx, sigmax);
  }
}
