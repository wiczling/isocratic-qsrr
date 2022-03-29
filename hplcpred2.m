function  samples = hplcpred2(samples,datastruct)

nSamples= size(samples.logkwHat,1);

nObs=datastruct.nObs;

samples.logkCond =zeros(nSamples, nObs);

for z = 1:nSamples
nu = samples.nu(z,:);
sigmaHat =samples.sigmaHat(z,:);
logkHat =samples.logkHat(z,:);
samples.logkCond(z,:)=logkHat+sigmaHat.*trnd(nu,size(logkHat));
end