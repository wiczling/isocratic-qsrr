%% The data and codes for 
% Title: 
% Authors: Pawe³ Wiczling, Agnieszka Kamedulska, £ukasz Kubik
% Adress: Department of Biopharmaceutics and Pharmacodynamics, Medical University of Gdañsk, Gen. J. Hallera 107, 80-416 Gdañsk, Poland
% Data: 30/10/2020
% Version 1.0
%% Load data
clc;
clear all;
data = readtable('Data\database_logk_1026.csv');
analyte_names = readtable('Data\database_logk_1026_analyte_names.csv');
functional_groups = readtable('Data\checkmol_functional_groups.csv');
functional_groups_names = readtable('Data\checkmol_functional_group_names.csv');

% combine nr of caroboxylic acid and carboxyalic acid salt functional groups
functional_groups{:,76}=functional_groups{:,76}+functional_groups{:,77};       
functional_groups{functional_groups{:,202}>8.1,202} = 8; % heterocyclic compounds with more than 8 heterocycles are treated as if they have 8 (strychnine)

% exclude functional groups that repeat itself (some groups are nested)
idx_excluded = [1 2 3 6 27 28 37 47 48 51 55 61 62 67 73 74 75 77 80 91 99 109 116 117 121 125 129 142 153 154 160 161 168 173 178 181 182 186 187 191 196];
writetable(functional_groups_names(idx_excluded,:),'Tables/functional_groups_excluded.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_excluded,:) = []; functional_groups(:,idx_excluded) = []; clear idx_excluded

% exclude functional groups not present on any analyte from the dataset
idx_not_present = find(sum(functional_groups{:,:})'==0);
writetable(functional_groups_names(idx_not_present,:),'Tables/functional_groups_not_present.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_not_present,:) = []; functional_groups(:,idx_not_present) = []; clear idx_not_present

%% Raw data
%Figure S1. Relationship between the logarithm of retention factor (log k)
%and acetonitrile content in the mobile phase. Lines connect measurements
%obtained for a particular analyte.

figure('Color',[1 1 1])
h1 = gscatter(data.concentration,data.logk,data.ID);
set(h1,'linestyle', '-')
xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
ylabel('logk')
legend off
box off
clear h1 

savefig('Figures/FigureS1_RawData.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS1_RawData.tif
%% Functional groups
% Figure 1. Functional groups identified by Checkmol. Figures show the
% number of analytes having at least one functional group of a given type.

[SortedSum,I] = sort(sum(functional_groups{:,:}>0.5));
figure('Color',[1 1 1])
subplot(1,2,1)
plot(1:1:50,SortedSum([1:1:50]),'-o')
xlabel('Functional group')
ylabel('                                                                      Number of analytes having at least one functional group of a given type')
view(90,90)
set(gca,'Xtick',[1:1:50],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([1:1:50]),2})
set(gca,'Yscale','lin','FontSize',8)
subplot(1,2,2)
plot(51:1:100,SortedSum([51:1:100]),'-o')
view(90,90)
set(gca,'Xtick',[51:1:100],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([51:1:100]),2})
set(gca,'Yscale','log','FontSize',8)
clear I SortedSum 

savefig('Figures/Figure1_FunctionalGroups.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure1_FunctionalGroups.tif
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5./2.1 18./2.1])
print -dtiff -r300 Figures/Figure1.tif

%% load logP and pKa
ACD_pKas = readtable([pwd '\Data\' 'ACD_pKas.csv']);
ACD_logP = readtable([pwd '\Data\' 'ACD_logP.csv']);
%% pKas
cov.pKaslit = ACD_pKas{:,3:5};
cov.pKasliterror = ACD_pKas{:,17:19};
cov.chargesA = abs(ACD_pKas{:,9:12});
cov.chargesB = abs(ACD_pKas{:,13:16});
cov.charges = cov.chargesA+cov.chargesB;  % absolute charge
cov.groupsA = diff(cov.chargesA,1,2);     % acidic group
cov.groupsB = -diff(cov.chargesB,1,2);    % basic group
cov.R = sum(cov.pKaslit<14,2);            % number of dissociation steps
cov.logP = ACD_logP.logP;                  % logP
cov.maxR = max(cov.R);
%% Initialize variables and parameters
clc
nObs = length(data.ID);
nAnalytes = length(unique(data.ID));
[~,~,j]=unique(data.ID,'first');

datastruct = struct(...
    'nObs',nObs, ...
	'nAnalytes', nAnalytes, ...
    'logPobs',cov.logP,...
    'K', size(functional_groups,2),...
    'nrfungroups',functional_groups{:,:},...
	'analyte',j,...
    'maxR',cov.maxR,...
    'R',cov.R,...
    'pKaslit',cov.pKaslit, ...
    'pKasliterror',cov.pKasliterror, ...
    'groupsA',cov.groupsA, ...
    'groupsB',cov.groupsB, ...
    'chargesA',cov.chargesA,...
    'chargesB',cov.chargesB,...
	'fi',data.concentration,...
    'run_estimation', 1, ...
	'logkObs', data.logk);
clear i1 j

%% Initialize
clear init0
% Initialize the values for each variable in each chain
for i=1:4
    S.logkwHat  =  normrnd(2.2,2,1);
	S.S1aHat     = normrnd(5,1,1) ;
    S.dlogkwHat = normrnd([-1 -1],0.125,1,2) ;
    S.dSaHat   = normrnd([0 0],0.5,1,2) ;
    S.S2aHat    = lognrnd(log(2),0.05,1,1) ; 
    S.beta  = normrnd([0.75 0.5],0.125,1,2) ;
    S.alphaAHt = normrnd([2 -1],0.2,1,2) ;
    S.sigma   = lognrnd(log(0.2),0.2,1, datastruct.nAnalytes);
    S.msigma  = lognrnd(log(0.2),0.2,1, 1);
    S.ssigma  = lognrnd(log(0.5),0.2,1, 1);
    S.omega = [1 1] .* exp(normrnd(0, 0.5, 1, 2));
    S.rho1 = [1 0.75
             0.75 1];
    S.kappa = [0.25 0.25] .* exp(normrnd(0, 0.2, 1, 2));
    S.tau   = [0.5] .* exp(normrnd(0, 0.2, 1, 1));
    S.pilogkw = zeros(1,datastruct.K);
    S.piS1a = zeros(1,datastruct.K);
    S.sdpi = [0.1 0.1] .* exp(normrnd(0, 0.1, 1, 2));
    S.param =  [2+1.*datastruct.logPobs 5*ones(datastruct.nAnalytes,1)+0.5.*datastruct.logPobs]; 
    S.etadlogkwA = zeros(datastruct.nAnalytes,datastruct.maxR+1);
    S.etadlogkwB = zeros(datastruct.nAnalytes,datastruct.maxR+1);
    S.etadSaA = zeros(datastruct.nAnalytes,datastruct.maxR+1);
    S.etadSaB = zeros(datastruct.nAnalytes,datastruct.maxR+1);
    S.pKaw = datastruct.pKaslit;
    S.alpha =zeros(datastruct.nAnalytes,datastruct.maxR);
    S.nu = gamrnd(2,1./0.1);
    init0(i) = S;
end
clear S i i1 j kaHat kwHat nAnalytes nObs fi nExp

%% 
fprintf( 'Running Stan...\n' );
fito= stan('file','hplc-qsrr-izo.stan','data', datastruct,'method','optimize', ...
              'working_dir','Tmpstan','verbose', logical(1),'init',init0(1),'iter',1000,'warmup',1000, ...  
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fito.block()
save('fito.mat', 'fito','-v7.3')
%% Use Stan.  
% Posterior 
datastruct.run_estimation=1;
fit= stan('file','hplc-qsrr-izo.stan','data', datastruct, 'verbose', logical(1), ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fit.block(); 
clc

%% Use Stan.  
% Prior. Better to generate them withouth running MCMC.
datastruct.run_estimation=0;
fitp= stan('file','hplc-qsrr-izo.stan','data', datastruct, 'verbose', logical(1), ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');      
%% Summary of model parameters. Save to file
% Table 1. Summary of the MCMC Simulations of the Marginal Posterior
% Distributions of Population-Level Model Parameters. Mean Denotes Sample
% Mean, MCSE Denotes Monte Carlo Standard Error, StdDev Denotes Sample
% Standard Deviation, 5%, 50%, 95% Denote Corresponding Quantiles, N_Eff
% Denotes Effective Sample Size, R_Hat Denotes a Measure of Chain
% Equilibrium, must be within 0.05 of 1.0.

diary hplc-qsrr-izo.txt
fit.print();
diary off
save hplc-qsrr-izo '-v7.3'
%% load saved data
load hplc-qsrr-izo
%% get samples
% samplesp = fitp.extract;
samples  = fit.extract;
samples = hplcpred2(samples,datastruct);


% %% Posterior summary (used later as priors)
% pmpilogkw = mean(samples.pilogkw)';
% pspilogkw = std(samples.pilogkw)';
% pmpidlogk = mean(samples.pidlogk)';
% pspidlogk = std(samples.pidlogk)';
% pmpilogS2 = mean(samples.pilogS2)';
% pspilogS2 = std(samples.pilogS2)';
% 
% for i=1:10
%     figure
%     subplot(3,1,1)
%     hold on
%     histogram(samples.pilogkw(:,i), 'Normalization','pdf'); 
%     plot(0:0.01:3,normpdf(-0:0.01:3,pmpilogkw(i),pspilogkw(i)))
%     subplot(3,1,2)
%     hold on
%     histogram(samples.pidlogk(:,i), 'Normalization','pdf'); 
%     plot(-1:0.01:3,normpdf(-1:0.01:3,pmpidlogk(i),pspidlogk(i)))
%     subplot(3,1,3)
%     hold on
%     histogram(samples.pilogS2(:,i), 'Normalization','pdf'); 
%     plot(-1:0.01:1,normpdf(-1:0.01:1,pmpilogS2(i),pspilogS2(i)))
% end
% 
% Posterior_summary = table(pmpilogkw,pspilogkw,pmpidlogk,pspidlogk,pmpilogS2,pspilogS2);
% 
% writetable(Posterior_summary,'Tables/Posterior_summary.csv','Delimiter',',','QuoteStrings',false)
% 
% clear pmpilogkw pspilogkw pmpilogka pspilogka pmpilogS2 pspilogS2
%% Prior predictive check (not used)
 % Visual predictive check:
 figure('Color', [1 1 1]);
 VPC(samplesp.logkCond', datastruct.logkObs, datastruct.fi, 1)
 title('Prior predicitve check')
% 
% Individaul and population predictions:
logkCond_p = prctile((samplesp.logkCond),[5 50 95],1);
logkPred_p = prctile((samplesp.logkCond),[5 50 95],1);
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:1026,10,'Replace',false);

for i = 1:10
    subplot(5, 2, i)
	hold on
	plot(datastruct.fi(datastruct.analyte==k(i)), datastruct.logkObs(datastruct.analyte==k(i)),  'k.')
	plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(1,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(2,datastruct.analyte==k(i)),  'k-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(3,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(1,datastruct.analyte==k(i)),  'r:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(2,datastruct.analyte==k(i)),  'r-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(3,datastruct.analyte==k(i)),  'r:')
    set(gca,'XTick',0:0.2:1)
    xlim([0 1])
	ylim([-6 6])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
         xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/PriorPredictionsModel.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/PriorPredictionsModel.tif

clear i logkPred_p logkCond_p
%% Goodness of Fit Plots, GOF
% Figure S7. Goodness-of-fit plots. The observed vs. the mean
% population-predicted retention factors (i.e., the a posteriori means of
% predictive distributions corresponding to the future observations of a
% new analyte) and the observed vs the mean individual-predicted retention
% times (i.e., the a posteriori mean of a predictive distribution
% conditioned on the observed data from the same analyte).
logkPred_mean  = mean(samples.logkCond);
logkCond_mean  = mean(samples.logkCond);
figure('Color', [1 1 1]);
subplot(2,1,1)
hold on
plot(logkPred_mean,datastruct.logkObs','.')
xlabel('Population predicted logk')
ylabel('Observed logk')
plot(xlim,xlim,':')
subplot(2,1,2)
hold on
plot(logkCond_mean,datastruct.logkObs','.')
plot(xlim,xlim,':')
xlabel('Individual Predicted logk')
ylabel('Observed logk')
clear logkPred_mean logkCond_mean

savefig('Figures/FigureS8_GOF.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS8_GOF.tif
%% Trace plots (not used)
samples_np = fit.extract('permuted',false);
Param = 'logkwHat';

[~,m]=size(samples_np(1).(Param));

for i=1:min(m,10)
figure('Color', [1 1 1]);
for z=1:4
hold on
plot(samples_np(z).(Param)(:,i),'-');
xlabel('Iteration')
end
ylabel([Param '(:,' num2str(i) ')'],'fontsize',12);
end

clear samples_np Param z i n m
%% Prior Posterior comparisons (not used):
Param = 'pilogkw';
[~,m]=size(samplesp.(Param));

for i=1:min(m,10)
figure('Color', [1 1 1]);
hold on
[f,xi] = ksdensity((samplesp.(Param)(:,i))); 
plot(xi,f);
[f,xi] = ksdensity((samples.(Param)(:,i)));  
plot(xi,f);
legend('Prior','Posterior')
xlabel([Param '(:,' num2str(i) ')'],'fontsize',12);
ylabel('Probability density estimate')
end

clear Param z i n m f xi 
%% Individaul and population predictions:
% Figure 3.  Individual and population predictions represented as posterior
% medians (lines) and 5th-95th percentiles (dotted lines) for a random set
% of 10 analytes. Observed retention factors are shown as dots. Black
% corresponds to future observations on the same analyte, and red
% corresponds to future observations of a new analyte.

logkCond_p = prctile((samples.logkCond),[5 50 95],1);
logkPred_p = prctile((samples.logkCond),[5 50 95],1);
figure('Color', [1 1 1]);
% rng(3333)
 k = datasample(1:1026,10,'Replace',false);

for i = 1:10
    subplot(5, 2, i)
	hold on
	plot(datastruct.fi(datastruct.analyte==k(i)), datastruct.logkObs(datastruct.analyte==k(i)),  'k.')
	plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(1,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(2,datastruct.analyte==k(i)),  'k-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(3,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(1,datastruct.analyte==k(i)),  'r:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(2,datastruct.analyte==k(i)),  'r-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(3,datastruct.analyte==k(i)),  'r:')
    set(gca,'XTick',0:0.2:1)
    xlim([0 1])
	ylim([-2.5 4])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
         xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/FigureS3_PredictionsModel.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS3_PredictionsModel.tif

clear i logkPred_p logkCond_p

%% Posterior predictive checks (not used)
figure('Color', [1 1 1]);
subset = datastruct.analyte>-inf;
VPC(samples.logkCond(:,subset)', datastruct.logkObs(subset), datastruct.fi(subset),0)
clear subset
%% Calculate individual (analyte-specific) parameters

idata.logP =cov.logP;
  for j = 1:1:datastruct.nAnalytes
      samples.etap(:,j,1) = samples.param(:,j, 1) - samples.miu(:,j,1); % 
      samples.etap(:,j,2) = samples.param(:,j, 2) - samples.miu(:,j,2); % 
  end
  
idata.etap = squeeze(mean(samples.etap));    % logkw, S1
idata.param = squeeze(mean(samples.param));  % logkw, S1

idata.logkwx = squeeze(mean(samples.logkwx));
idata.logkax = squeeze(mean(samples.logkwx-samples.S1ax));


idata.pKaw = squeeze(mean(samples.pKaw));
idata.alphaa = squeeze(mean(samples.alphaa));

idata.pKaa=squeeze(mean(samples.pKaw+samples.alphaa));

  for j = 1:1:datastruct.nAnalytes
 samples.dlogkw(:,j,:) =  (squeeze(samples.dlogkwHat(:,1))*ones(1,4) + (squeeze(samples.kappa(:,1))*ones(1,4)).*(squeeze(samples.etadlogkwA(:,j,:)))) .* datastruct.chargesA(j,:) ...
                        + (squeeze(samples.dlogkwHat(:,2))*ones(1,4) + (squeeze(samples.kappa(:,1))*ones(1,4)).*(squeeze(samples.etadlogkwB(:,j,:)))) .* datastruct.chargesB(j,:);
 samples.dS1a(:,j,:) =    (squeeze(samples.dSaHat(:,1))*ones(1,4)    + (squeeze(samples.kappa(:,2))*ones(1,4)).*squeeze(samples.etadSaA(:,j,:))  )    .* datastruct.chargesA(j,:) ...
                        + (squeeze(samples.dSaHat(:,2))*ones(1,4)    + (squeeze(samples.kappa(:,2))*ones(1,4)).*squeeze(samples.etadSaB(:,j,:)) )     .* datastruct.chargesB(j,:) ;
  end

 samples.dlogka = samples.dlogkw - samples.dS1a;

 idata.dlogkw = squeeze(mean(samples.dlogkw))./cov.charges;
 idata.dlogka = squeeze(mean(samples.dlogka))./cov.charges;
 idata.dS1a   = squeeze(mean(samples.dS1a))  ./cov.charges;
 
 idata.chargesAB = cov.chargesB - cov.chargesA;   % {-2,-1,0,1,2}
 idata.groupsAB  = cov.groupsB  - cov.groupsA;    % {-1 for Acids, 1 for Bases}

 idata.isdiss = 0.*cov.charges;
 for j = 1:1:datastruct.nAnalytes
    idata.isdiss(j,1:cov.R(j)+1)=1;
 end

clear j

%% Individual Parameters - Neutral Form
figure('Color', [1 1 1]);
xynames = {'logkwN_{i}','S1aN_{i}','logP_i'};
gplotmatrix([idata.param(:,1) idata.param(:,2) idata.logP],[],0*idata.logP,'kk',[],[],'on','stairs',xynames,xynames)

 h=get(gcf,'children');
 set(h(1),'Visible','off')
 savefig('Figures/IndividualParametersNeutralForm.fig')
 set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
 print -dtiff -r300 Figures/IndividualParametersNeutralForm.tif
clear h xynames
%% Effect of dissociation
figure('Color', [1 1 1]);
xynames = {'dlogkw_{r,i}','dS1a_{r,i}'};
ktore = idata.isdiss(:)~=0;
gplotmatrix([idata.dlogkw(ktore) idata.dS1a(ktore)],[],idata.chargesAB(ktore),'rrrybbb',[],[],'on','stairs',xynames,xynames)

h=get(gcf,'children');
set(h(1),'Visible','off')
savefig('Figures/IndParamEffectDiss.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/IndParamEffectDiss.tif
clear h xynames ktore

%% Influence of functional groups
% Figure 2. Graphical display of the marginal posterior distributions for
% the effects of each functional group on logkw, logka, and logS2.
figure('Color', [1 1 1]);
subplot(1,4,2)
hold on
boxplot_pwhisker(samples.pilogkw(:,:),{'Labels',functional_groups_names{:,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2139    0.1100    0.2138    0.8150])
ylabel('\pi_{logk_w}','FontSize',8)
subplot(1,4,3)
hold on
boxplot_pwhisker(samples.piS1a(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1.5])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.4854    0.1100    0.2178    0.8150])
ylabel('\pi_{S1_a}','FontSize',8)
subplot(1,4,4)
hold on
% boxplot_pwhisker(samples.pilogS2(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-0.5 0.5])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.7334    0.1100    0.1708    0.8150])
ylabel('\pi_{logS_2}','FontSize',8)

savefig('Figures/Figure2_FunctionalGroupEffects.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2_FunctionalGroupEffects.tif

%% Influence of functional groups (pilogkw)
% Figure S2A. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on pilogkw.
figure('Color', [1 1 1]);
[~,I] = sort(mean(samples.pilogkw(:,:)),'descend');
hold on
boxplot_pwhisker(samples.pilogkw(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{logk_w}','FontSize',8)

savefig('Figures/Figure2A_FunctionalGroupEffects_pilogkw.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2A_FunctionalGroupEffects_pilogkw.tif
%% Influence of functional groups (piS1a)
% Figure S2B. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on piS1a.
figure('Color', [1 1 1]);

[~,I] = sort(mean(samples.piS1a(:,:)),'descend');

hold on
boxplot_pwhisker(samples.piS1a(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 2])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{piS1a}','FontSize',8)

savefig('Figures/Figure2B_FunctionalGroupEffects_piS1a.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2B_FunctionalGroupEffects_piS1a.tif

%% Influence of functional groups (pilogka)
% Figure S2C. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on pilogS2.

figure('Color', [1 1 1]);
[~,I] = sort(mean(samples.pilogkw(:,:)-samples.piS1a(:,:)),'descend');

hold on
boxplot_pwhisker(samples.pilogkw(:,I)-samples.piS1a(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-2 2])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{logka}','FontSize',8)

savefig('Figures/Figure2C_FunctionalGroupEffects_pilogka.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2C_FunctionalGroupEffects_pilogka.tif

%%
ver
% ----------------------------------------------------------------------------------------------------
% MATLAB Version: 9.2.0.556344 (R2017a)
% MATLAB License Number: 261217
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 19041)
% Java Version: Java 1.7.0_60-b19 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% ----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.2         (R2017a)
% Bioinformatics Toolbox                                Version 4.8         (R2017a)
% Curve Fitting Toolbox                                 Version 3.5.5       (R2017a)
% Global Optimization Toolbox                           Version 3.4.2       (R2017a)
% MATLAB Compiler                                       Version 6.4         (R2017a)
% MATLAB Compiler SDK                                   Version 6.3.1       (R2017a)
% Optimization Toolbox                                  Version 7.6         (R2017a)
% Parallel Computing Toolbox                            Version 6.10        (R2017a)
% SimBiology                                            Version 5.6         (R2017a)
% Statistics and Machine Learning Toolbox               Version 11.1        (R2017a)
% Symbolic Math Toolbox                                 Version 7.2         (R2017a)

%% Licenses 
%1) Code &copy; 2020, Pawe³ Wiczling, licensed under BSD-3.
%2) Text &copy; 2020, Pawe³ Wiczling, licensed under CC-BY-NC 4.0.