functions {

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
vector lower_tri(matrix mat) {

int d = rows(mat);
int lower_tri_d = d * (d - 1) / 2;
vector[lower_tri_d] lower;
int count = 1;
for(r in 2:d) {
for(c in 1:(r - 1)) {
lower[count] = mat[r,c];
count += 1;
}
}
return(lower); 
}

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
real lkj_corr_point_lower_tri_lpdf(matrix rho, real point_mu_lower, real point_scale_lower) {

real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(rho[2,1] | point_mu_lower, point_scale_lower);
return(lpdf);
}

real funlogki(vector logkwx, 
           vector S1x, 
           vector pKaw, 
           vector alpha,
           real S2,
           int nDiss,
           real fi) {

real logki;
vector[4] S1;
vector[4] logkix;
vector[3] pHpKa;


logkix   = logkwx-S1x*fi/(1+S2*fi);

pHpKa    = 2.6-(pKaw+alpha*fi);

if (nDiss==0) {
  logki = logkix[1]; 
}

else if (nDiss==1){
  logki = logkix[1] +
          log1p_exp(log(10)*(pHpKa[1]+logkix[2]-logkix[1]))/log(10)-
          log1p_exp(log(10)*(pHpKa[1]))/log(10);
}

else if (nDiss==2){
  logki = logkix[1] +
          log1p_exp(log(10)*(pHpKa[1]+logkix[2]-logkix[1]) + 
          log1p_exp(log(10)*(pHpKa[2]+logkix[3]-logkix[2])))/log(10)-
          log1p_exp(log(10)*(pHpKa[1]) + 
          log1p_exp(log(10)*(pHpKa[2])))/log(10);
}
else if (nDiss==3){
  logki = logkwx[1] +
          log1p_exp(log(10)*(pHpKa[1]+logkix[2]-logkix[1]) + 
          log1p_exp(log(10)*(pHpKa[2]+logkix[3]-logkix[2]) + 
          log1p_exp(log(10)*(pHpKa[3]+logkix[4]-logkix[3]))))/log(10)-
          log1p_exp(log(10)*(pHpKa[1]) + 
          log1p_exp(log(10)*(pHpKa[2]) + 
          log1p_exp(log(10)*(pHpKa[3]))))/log(10);
}
else if (nDiss==4){
  logki = logkwx[1] +
          log1p_exp(log(10)*(pHpKa[1]+logkix[2]-logkix[1]) + 
          log1p_exp(log(10)*(pHpKa[2]+logkix[3]-logkix[2]) + 
          log1p_exp(log(10)*(pHpKa[3]+logkix[4]-logkix[3]) + 
          log1p_exp(log(10)*(pHpKa[4]+logkix[5]-logkix[4])))))/log(10)-
          log1p_exp(log(10)*(pHpKa[1]) + 
          log1p_exp(log(10)*(pHpKa[2]) + 
          log1p_exp(log(10)*(pHpKa[3]) + 
          log1p_exp(log(10)*(pHpKa[4])))))/log(10);

}
else if (nDiss==5){
  logki = logkwx[1] +
          log1p_exp(log(10)*(pHpKa[1]+logkix[2]-logkix[1]) + 
          log1p_exp(log(10)*(pHpKa[2]+logkix[3]-logkix[2]) + 
          log1p_exp(log(10)*(pHpKa[3]+logkix[4]-logkix[3]) + 
          log1p_exp(log(10)*(pHpKa[4]+logkix[5]-logkix[4]) + 
          log1p_exp(log(10)*(pHpKa[5]+logkix[6]-logkix[5]))))))/log(10)-
          log1p_exp(log(10)*(pHpKa[1]) + 
          log1p_exp(log(10)*(pHpKa[2]) + 
          log1p_exp(log(10)*(pHpKa[3]) + 
          log1p_exp(log(10)*(pHpKa[4]) + 
          log1p_exp(log(10)*(pHpKa[5]))))))/log(10);
}

return logki;
}

}

data{
int nAnalytes;	           // number of analytes
int nObs;		           // number of observations
int analyte[nObs];	       // analyte indexes
real fi[nObs]; 
vector[nAnalytes] logPobs; 
int<lower=0,upper=5> maxR;
int<lower=0,upper=5> R[nAnalytes];
ordered[maxR] pKaslit[nAnalytes];
vector[maxR] pKasliterror[nAnalytes];
vector[maxR] groupsA[nAnalytes];
vector[maxR] groupsB[nAnalytes];
vector[maxR+1] chargesA[nAnalytes];
vector[maxR+1] chargesB[nAnalytes];
int<lower=0> K;                      //  number of predictors (functional groups)
matrix[nAnalytes, K] nrfungroups;    // predictor matrix (functional groups)   
vector[nObs] logkObs; // observed logkObs
int<lower = 0, upper = 1> run_estimation; // 0 for prior predictive, 1 for estimation 
}

transformed data {

}

parameters{
real logkwHat;	       // typical value of logkw [N]
real S1aHat;           // typical value of S1a [N]
real dlogkwHat[2];     // typical value of dlogkw [A,B] 
real dSaHat[2];        // typical value of dlogka [A,B] 
real<lower = 0> S2aHat; // typical value of S2a
vector[2] beta;         // effects of logP 
vector[2] alphaHat;  // changes of pKa with org. mod 
vector<lower = 0.01>[2] omega;   // between analyte variabilities (neutral forms)
corr_matrix[2] rho1;	                        // correlation matrix	 
vector<lower = 0.01>[2] kappa;    // between analyte variabilities (diss. forms)
real<lower = 0.01> tau;     // between analyte variabilities for acids pKa

vector[K] pilogkw;  // regression coefficient for logkw
vector[K] piS1a;  // regression coefficient for S1a

vector<lower = 0.01>[2] sdpi;     // between analyte variabilities for acids pKa

// residual variability
real<lower = 0.01> msigma; // mean
real<lower = 0.01> ssigma; // scale
real<lower = 1.1> nu;
// individual values of chromatographic parameters
vector[2] param[nAnalytes]; 	
matrix[nAnalytes,maxR+1] etadlogkwA;
matrix[nAnalytes,maxR+1] etadlogkwB;
matrix[nAnalytes,maxR+1] etadSaA;
matrix[nAnalytes,maxR+1] etadSaB;
matrix[nAnalytes,maxR] alpha; 
vector[maxR] pKaw[nAnalytes];

// and residuals
vector<lower = 0.01, upper = 4>[nAnalytes] sigma;
}

transformed parameters{
vector[maxR+1] logkwx[nAnalytes];
vector[maxR+1] S1mx[nAnalytes];
vector[maxR+1] S1ax[nAnalytes];
vector[maxR] alphaa[nAnalytes];

vector[2] miu[nAnalytes];	
cov_matrix[2] Omega; // variance-covariance matrix

vector[nObs] logkHat;
vector[nObs] sigmaHat;

Omega = quad_form_diag(rho1, omega);	// diag_matrix(omega) * rho * diag_matrix(omega)

for(i in 1:nAnalytes){
    miu[i,1]  = logkwHat + beta[1] * (logPobs[i]-2.2) + nrfungroups[i,1:K] * pilogkw;
    miu[i,2]  = S1aHat   + beta[2] * (logPobs[i]-2.2) + nrfungroups[i,1:K] * piS1a;
}

for(i in 1:nAnalytes){
for(r in 1:maxR+1){
logkwx[i,r] = param[i, 1] +
            (dlogkwHat[1]+kappa[1]*etadlogkwA[i,r])*chargesA[i,r] +
            (dlogkwHat[2]+kappa[1]*etadlogkwB[i,r])*chargesB[i,r];
S1ax[i,r] = (param[i, 2] + 
            (dSaHat[1]+kappa[2]*etadSaA[i,r])*chargesA[i,r] +
            (dSaHat[2]+kappa[2]*etadSaB[i,r])*chargesB[i,r])*(1+S2aHat);
}}

for(i in 1:nAnalytes){
for(r in 1:maxR){
alphaa[i,r] = (alphaHat[1]+tau*alpha[i,r]) * groupsA[i,r] + (alphaHat[2]+tau*alpha[i,r])*groupsB[i,r];
}}


for(i in 1:nObs){

 logkHat[i] = funlogki(logkwx[analyte[i],], 
               S1ax[analyte[i],],  
               pKaw[analyte[i],], 
               alphaa[analyte[i],],
               S2aHat,
               R[analyte[i]],
               fi[i]);
 
  sigmaHat[i] = sigma[analyte[i]];
}


}
model{
logkwHat  ~ normal(2.2, 2);
S1aHat    ~ normal(5, 1);
dlogkwHat ~ normal(-1,0.125);
dSaHat    ~ normal(0,0.5);
S2aHat    ~ lognormal(0.69,0.125);
alphaHat[1] ~ normal(2,0.25);
alphaHat[2] ~ normal(-1,0.25);
beta[1] ~ normal(1,0.125);
beta[2] ~ normal(0.5,0.5);
omega     ~ normal(0,2);
rho1       ~ lkj_corr_point_lower_tri(0.75, 0.125);
kappa      ~ normal(0,0.5);

pilogkw ~ normal(0,sdpi[1]);
piS1a   ~ normal(0,sdpi[2]);
sdpi ~ normal(0,0.1);
tau ~ normal(0,0.5);

sigma  ~ lognormal(log(msigma),ssigma); 
msigma ~ normal(0,1);
ssigma ~ normal(0,1);

nu ~ gamma(2,0.1);

for(i in  1:nAnalytes){
param[i] ~ multi_normal(miu[i],Omega);
}

to_vector(etadlogkwA) ~ std_normal();
to_vector(etadlogkwB) ~ std_normal();
to_vector(etadSaA) ~ std_normal();
to_vector(etadSaB) ~ std_normal();

to_vector(alpha) ~ std_normal();

for (i in 1:nAnalytes) pKaw[i] ~ normal(pKaslit[i],pKasliterror[i]);

  if(run_estimation==1){
  logkObs ~ student_t(nu,logkHat,sigmaHat);
  }
}

generated quantities{
}