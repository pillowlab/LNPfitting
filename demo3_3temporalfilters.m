% demo3_3temporalfilters.m
%
% Tutorial script illustrating maximum likelihood / maximally informative
% dimensions (MID) estimation for LNP model with THREE temporal filters
% using both CBF and RBF nonlinearities.
%
% - also computes iSTAC estimate of three-filter model, for comparison.
% - Compares filter estimates using basis reconstructions of true filters
% - Computes test log-likelihood in units of bits / spike (equal to
%       information gain over a homogeneous Poisson model)
% - plots rate predictions of three different models (iSTC, ml-cbf, ml-rbf)
%       on test data. 

% initialize paths
initpaths;

datasetnum = 1;  % select: 1 (white noise) or 2 (correlated)
trainfrac = .8; % fraction of data to use for training (remainder is "test data")

% Load data divided into training and test sets
[Stim_tr,sps_tr,Stim_tst,sps_tst,RefreshRate,filts_true] = loadSimDataset(datasetnum,trainfrac);

slen_tr = size(Stim_tr,1);   % length of training stimulus / spike train
slen_tst = size(Stim_tst,1); % length of test stimulus / spike train
nsp_tr = sum(sps_tr);   % number of spikes in training set
nsp_tst = sum(sps_tst); % number of spikes in test set


%% == 2. Compute iSTAC estimator (for comparison sake)

nkt = 30; % number of time bins to use for filter 
% This is important: normally would want to vary this to find optimal filter length

% Compute STA and STC
[sta,stc,rawmu,rawcov] = simpleSTC(Stim_tr,sps_tr,nkt);  % compute STA and STC

% Compute iSTAC estimator
nistacFilts = 3; % number of iSTAC filters to compute
fprintf('\nComputing iSTAC estimate\n');
[istacFilts,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nistacFilts); % find iSTAC filters

% Fit iSTAC model nonlinearity using filters 1 and 2
pp_istac = fitNlin_expquad_ML(Stim_tr,sps_tr,istacFilts,RefreshRate); % LNP model struct
% Compute training and test log-likelihood
LListac_tr = logli_LNP(pp_istac,Stim_tr,sps_tr); % training log-likelihood
[LListac_tst,rate_istac] = logli_LNP(pp_istac,Stim_tst,sps_tst); % test log-likelihood

% -----  Plot filters ------------
subplot(231); 
tt = -nkt+1:0;

% Compute true filters reconstructed in basis of iSTAC estimates
filtsHat_istac = istacFilts*(istacFilts\filts_true);

plot(tt,filts_true,'k--',tt,filtsHat_istac,'linewidth',2); 
axis tight; title('Filter reconstructions: istac estimates');
xlabel('time before spike (bins)');


%% == 3. Set up temporal basis for stimulus filters (for ML / MID estimators)

pp0 = makeFittingStruct_LNP(istacFilts(:,1),RefreshRate); % initialize param struct

% == Set up temporal basis for representing filters (see demo 1 more detail)  ====
% (try changing these params until basis can accurately represent all iSTAC axes).
ktbasprs.neye = 0; % number of "identity"-like basis vectors
ktbasprs.ncos = 12; % number of raised cosine basis vectors
ktbasprs.kpeaks = [1 2*nkt/3]; % location of 1st and last basis vector bump
ktbasprs.b = 1.5; % determines how nonlinearly to stretch basis (higher => more linear)
[ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, nkt); % make basis
filtprs_basis = (ktbas'*ktbas)\(ktbas'*sta);  % filter represented in new basis
sta_basis = ktbas*filtprs_basis;

% Insert filter basis into fitting struct
pp0.k = sta_basis; % insert sta filter
pp0.kt = filtprs_basis; % filter coefficients (in temporal basis)
pp0.ktbas = ktbas; % temporal basis
pp0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

%% == 4. ML/MID 1:  ML estimator for LNP with CBF (cylindrical basis func) nonlinearity

nfilts = 3; % number of filters to recover

% Set parameters for cylindrical basis funcs (CBFs) and initialize fit
fstructCBF.nfuncs = 3; % number of basis functions for nonlinearity
fstructCBF.epprob = [.01, 0.99]; % cumulative probability outside outermost basis function peaks
fstructCBF.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function
fstructCBF.nlfuntype = 'cbf';

% Fit the model (iteratively adding one filter at a time)
[ppcbf,negLcbf] = fitLNP_multifiltsCBF_ML(pp0,Stim_tr,sps_tr,nfilts,fstructCBF,istacFilts);

% Compute test log-likelihood
LLcbf_tr = logli_LNP(ppcbf,Stim_tr,sps_tr); % training log-likelihood
[LLcbf_tst,rate_cbf] = logli_LNP(ppcbf,Stim_tst,sps_tst); % test log-likelihood

% Compute true filters reconstructed in basis of iSTAC estimates
filts_cbf = squeeze(ppcbf.k);  % filter estimates
filtsHat_cbf = filts_cbf*(filts_cbf\filts_true); % reconstructed true filts

% Plot reconstruction of true filters
subplot(232);
plot(tt,filts_true,'k--',tt,filtsHat_cbf,'linewidth',2); 
axis tight; title('Filter reconstructions: ML-cbf estimates');


%% == 5. ML / MID 2:  ML estimator for LNP with RBF (radial basis func) nonlinearity

% Set parameters for cylindrical basis funcs (CBFs) and initialize fit
fstructRBF.nfuncs = 3; % number of basis functions for nonlinearity
fstructRBF.epprob = [.01, 0.99]; % cumulative probability outside outermost basis function peaks
fstructRBF.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function
fstructRBF.nlfuntype = 'rbf';

% Fit the model (iteratively adding one filter at a time)
[pprbf,negLrbf] = fitLNP_multifiltsRBF_ML(pp0,Stim_tr,sps_tr,nfilts,fstructRBF,istacFilts);

% Compute test log-likelihood
LLrbf_tr = logli_LNP(pprbf,Stim_tr,sps_tr); % training log-likelihood
[LLrbf_tst,rate_rbf] = logli_LNP(pprbf,Stim_tst,sps_tst); % test log-likelihood

% Compute true filters reconstructed in basis of iSTAC estimates
filts_rbf = squeeze(pprbf.k);  % filter estimates
filtsHat_rbf = filts_rbf*(filts_rbf\filts_true); % reconstructed true filts

% Plot reconstruction of true filters
subplot(233);
plot(tt,filts_true,'k--',tt,filtsHat_rbf,'linewidth',2); 
axis tight; title('Filter reconstructions: ML-rbf estimates');


%% 6. ==== Compute training and test performance in bits/spike =====


% Let's put istac results on same footing by representing them in the same
% temporal basis

% Compute true filters reconstructed in basis of iSTAC estimates
istacFiltsBasis = ktbas*(ktbas\istacFilts);
filtsHat_istacBasis = istacFiltsBasis*(istacFiltsBasis\filts_true);

% redo first plot
subplot(231);
plot(tt,filts_true,'k--',tt,filtsHat_istacBasis,'linewidth',2); 
axis tight; title('Filter reconstructions: istac estimates (in kt basis)');
xlabel('time bin before spike');

% ==== report R2 error in reconstructing first two filters =====
kmse = sum((filts_true(:)-mean(filts_true(:)).^2)); % mse of these two filters around mean
ferr = @(k)(sum((filts_true(:)-k(:)).^2)); % error in optimal R2 reconstruction of ktrue
fRsq = @(k)(1-ferr(k)/kmse); % R-squared

Rsqvals = [fRsq(filtsHat_istacBasis),fRsq(filtsHat_cbf),fRsq(filtsHat_rbf)];

fprintf('\n=========== RESULTS ======================================\n');
fprintf('\nFilter R^2:\n------------\n');
fprintf('istac:%.2f cbf:%.2f rbf:%.2f\n',Rsqvals);

% Plot filter R^2 values
subplot(245);
axlabels = {'istac','cbf','rbf'};
bar(Rsqvals); ylabel('R^2'); title('filter R^2');
set(gca,'xticklabel',axlabels, 'ylim', [.9 1]);


% ====== Compute log-likelihood / single-spike information ==========

% Compute log-likelihood of constant rate (homogeneous Poisson) model
muspike_tr = nsp_tr/slen_tr;    % mean number of spikes / bin, training set
muspike_tst = nsp_tst/slen_tst; % mean number of spikes / bin, test set
LL0_tr =   nsp_tr*log(muspike_tr) - slen_tr*muspike_tr; % log-likelihood, training data
LL0_tst = nsp_tst*log(muspike_tst) - slen_tst*muspike_tst; % log-likelihood test data

% Functions to compute single-spike informations
f1 = @(x)((x-LL0_tr)/nsp_tr/log(2)); % compute training single-spike info 
f2 = @(x)((x-LL0_tst)/nsp_tst/log(2)); % compute test single-spike info
% (Divide by log 2 to get 'bits' instead of 'nats')

SSinfo_tr = [f1(LListac_tr), f1(LLcbf_tr) f1(LLrbf_tr)];
SSinfo_tst = [f2(LListac_tst),f2(LLcbf_tst) f2(LLrbf_tst)];

fprintf('\nSingle-spike information (bits/spike):\n');
fprintf('------------------------------------- \n');
fprintf('Train: istac: %.2f  cbf:%.2f  rbf:%.2f\n', SSinfo_tr);
fprintf('Test:  istac: %.2f  cbf:%.2f  rbf:%.2f\n', SSinfo_tst);

% Plot test single-spike information
subplot(246);
bar(SSinfo_tst); ylabel('bits / sp'); title('test single spike info');
set(gca,'xticklabel',axlabels, 'ylim', [0.75 1]);

% Plot some rate predictions for first 100 bins
subplot(2,4,7:8);
ii = 1:100;
plot(ii,rate_istac(ii),ii,rate_cbf(ii),ii,rate_rbf(ii), 'linewidth',2);
title('rate predictions on test data (3-filter models)');
legend('istac', 'ml-cbf','ml-rbf');
ylabel('rate (sp/s)'); xlabel('time bin');