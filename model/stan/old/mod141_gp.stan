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

  real lkj_corr_cholesky_point_lower_tri_lpdf(matrix rhoc, real point_mu_lower, real point_scale_lower) {
// works only for [2x2 matrix]
    matrix[2,2] rho = multiply_lower_tri_self_transpose(rhoc);
    real lpdf = lkj_corr_cholesky_lpdf(rhoc | 1);
    lpdf = + normal_lpdf(rho[2, 1] | point_mu_lower, point_scale_lower);
    return lpdf;
  }
  
matrix chol_kronecker_prod(matrix A, matrix B) {
  //https://discourse.mc-stan.org/t/stan-efficient/26329/8
int m = rows(A);
int p = rows(B);
int N = m * p;
matrix[N, N] C = rep_matrix(0., N, N);

for (i in 1:m) {
for (j in 1:i) {
int row_start = (i - 1) * p + 1;
int row_end = (i - 1) * p + p;
int col_start = (j - 1) * p + 1;
int col_end = (j - 1) * p + p;
C[row_start:row_end, col_start:col_end] = A[i, j] * B;
}
}
return C;
}

  matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

  real funlogki2(vector logkw, vector S1, real S2, vector pKaw, vector alpha, int R, real fi) {
    real lnki;
    real logki;
    vector[4] logkix;
    vector[3] pHmpKa;
    
    logkix = log(10) *(logkw - S1*(1+S2) * fi / (1 + S2 * fi));
    pHmpKa = log(10) *(2.66 - (pKaw + alpha * fi));
    
    if (R == 0) {
    lnki = logkix[1];
    } else if (R == 1) {
    lnki =  log_sum_exp([logkix[1],
                        logkix[2]+pHmpKa[1]])-
            log_sum_exp([0.0,
                        pHmpKa[1]]);
    } else if (R == 2) {
    lnki =  log_sum_exp([logkix[1],
                        logkix[2]+pHmpKa[1],
                        logkix[3]+pHmpKa[1]+pHmpKa[2]])-
            log_sum_exp([0.0,
                        pHmpKa[1],
                        pHmpKa[1]+pHmpKa[2]]);
    } else if (R == 3) {
    lnki =  log_sum_exp([logkix[1],
                        logkix[2]+pHmpKa[1],
                        logkix[3]+pHmpKa[1]+pHmpKa[2],
                        logkix[4]+pHmpKa[1]+pHmpKa[2]+pHmpKa[3]])-
            log_sum_exp([0.0,
                        pHmpKa[1],
                        pHmpKa[1]+pHmpKa[2],
                        pHmpKa[1]+pHmpKa[2]+pHmpKa[3]]);
             
        
    }
    
    logki = lnki/log(10);
    
    return logki;
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
  int<lower=0> nK;                       // number of predictors (functional groups)
  matrix[nAnalytes, nK] fgrp;     // predictor matrix (functional groups) 
  vector[nObs] logkobs;                 // observed retention factors 
  cov_matrix[nAnalytes] covK; //  
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  
  array[nAnalytes] vector[maxR] alphax;
  array[nAnalytes] vector[maxR] pKawx;
  
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
  
}
parameters {
  real logkwHat;            // typical logkw [Neutral]
  real S1Hat;               // effect of ACN on logkw [Neutral]
  real dlogkwHat;  // effect of dissociation on logkw [Acids, Bases]
  real dS1Hat;     // effect of dissociation on S1 [Acids, Bases] 
  real<lower=0> S2Hat;            // typical value of S2a (log10 scale)
  vector[2] beta;           // effect of logP on logkw and S1
  vector<lower=0>[2] omega; // sd of BAV [logkw,S1]
  corr_matrix[2] rho;       // correlation matrix [logkw vs. S1] 
  real<lower=0> sigma;      // typical sigma
  vector[nK] pilogkw;        // pi-logkw
  vector[nK] piS1;           // pi-S1
  vector<lower = 0.01>[2] sdpi;     // between analyte variabilities for fgrp
  matrix[nAnalytes,2] paramN;
}

transformed parameters {
  matrix[nAnalytes,2] miu;
  array[nAnalytes] vector[maxR + 1] logkwx;
  array[nAnalytes] vector[maxR + 1] S1x;
  vector[nObs] logkx;
  cov_matrix[2] Omega;
  matrix[nAnalytes*2,nAnalytes*2] OmegaK; 
  
  Omega = quad_form_diag(rho, omega); // diag_matrix(omega) * rho * diag_matrix(omega)
  OmegaK = chol_kronecker_prod(Omega,covK);
    
  for (i in 1 : nAnalytes) {
    miu[i, 1] = logkwHat + beta[1] * (logPobs[i]-2.5) + fgrp[i,1:nK] * pilogkw;
    miu[i, 2] = S1Hat    + beta[2] * (logPobs[i]-2.5) + fgrp[i,1:nK] * piS1;
  }
  
  for (d in 1 : nA0) {
     miu[idxA0[d], 1] += dlogkwHat;
     miu[idxA0[d], 2] += dS1Hat ;
   }
    for (d in 1 : nB0) {
     miu[idxB0[d], 1] += dlogkwHat ;
     miu[idxB0[d], 2] += dS1Hat ;
   }
    
  for (i in 1 : nAnalytes) { 
   logkwx[i, : ]  =  paramN[i,1]*[1,1,1,1]';
      S1x[i, : ]  =  paramN[i,2]*[1,1,1,1]';
  }

  for (d in 1 : nGroupsA) {
      logkwx[idxGroupsA[d,1],  (idxGroupsA[d,2]+1) : 4] += dlogkwHat ;
         S1x[idxGroupsA[d,1],  (idxGroupsA[d,2]+1) : 4] += dS1Hat ;
         }
  
  for (d in 1 : nGroupsB) {
      logkwx[idxGroupsB[d,1], 1:idxGroupsB[d,2]] += dlogkwHat ;
      S1x[idxGroupsB[d,1],    1:idxGroupsB[d,2]] += dS1Hat ;
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
 logkwHat ~ normal(4, 4);
  S1Hat ~ normal(4, 2);
  dlogkwHat ~ normal(-1, 0.25);
  dS1Hat ~ normal(0, 0.25);
  S2Hat ~ lognormal(0.693, 0.125);
  beta[{1}] ~ normal(0.7, 0.125);
  beta[{2}] ~ normal(0.5, 0.5);
  omega ~ normal(0, 1);
  
  pilogkw ~ normal(0,sdpi[1]);
  piS1    ~ normal(0,sdpi[2]);
  sdpi ~ normal(0,0.1);
  
  rho ~ lkj_corr_point_lower_tri(0.75, 0.125);
 
  to_vector(paramN) ~ multi_normal_cholesky(to_vector(miu), OmegaK);
  
  sigma ~ normal(0,0.05);
  
  if (run_estimation == 1) {
  logkobs ~ normal(logkx, sigma);
  }
}

generated quantities {

}

