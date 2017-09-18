% demo2_2temporalfilters.m
%
% Tutorial script illustrating maximum likelihood / maximally informative
% dimensions (MID) estimation for LNP model with TWO temporal filters

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
fprintf('\nComputing iSTAC estimate\n');
nistacFilts = 3; % number of filters to compute
[istacFilts,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nistacFilts); % find iSTAC filters

% Fit iSTAC model nonlinearity using filters 1 and 2
pp_istac12 = fitNlin_expquad_ML(Stim_tr,sps_tr,istacFilts(:,1:2),RefreshRate); % LNP model struct

% Compute training and test log-likelihood
LListac12_tr = logli_LNP(pp_istac12,Stim_tr,sps_tr); % training log-likelihood
LListac12_tst = logli_LNP(pp_istac12,Stim_tst,sps_tst); % test log-likelihood

% Fit iSTAC model nonlinearity using filters 2 and 3
pp_istac23 = fitNlin_expquad_ML(Stim_tr,sps_tr,istacFilts(:,2:3),RefreshRate); % LNP model struct

% Compute training and test log-likelihood
LListac23_tr = logli_LNP(pp_istac23,Stim_tr,sps_tr); % training log-likelihood
LListac23_tst = logli_LNP(pp_istac23,Stim_tst,sps_tst); % test log-likelihood

% -----  Visualize filters and 2D nonlinearity for both models ---------------------
[xgrd12,ygrd12,nlvals12] = compNlin_2D(pp_istac12.k,pp_istac12.nlfun,Stim_tr); % compute gridded 2D nonlinearity
[xgrd23,ygrd23,nlvals23] = compNlin_2D(pp_istac23.k,pp_istac23.nlfun,Stim_tr); % compute gridded 2D nonlinearity

subplot(241); 
tt = -nkt+1:0;
plot(tt,filts_true,'k--',tt,istacFilts(:,1:2),'linewidth',2); 
axis tight; title('istac Filters 1 & 2');

subplot(245);
mesh(xgrd12,ygrd12,nlvals12); 
zlm = [0 max(nlvals12(:))*1.01];
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity: filts 1 & 2')

subplot(242); 
plot(tt,filts_true,'k--',tt,istacFilts(:,2:3),'linewidth',2); 
axis tight; title('istac Filters 2 & 3');

subplot(246); mesh(xgrd23,ygrd23,nlvals23); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 2 axis');ylabel('filter 3 axis'); 
title('2D iSTAC nonlinearity: filts 2 & 3')

%% == 3. Set up temporal basis for stimulus filters (for ML / MID estimators)

pp0 = makeFittingStruct_LNP(istacFilts(:,1),RefreshRate); % initialize param struct

% == Set up temporal basis for representing filters (see demo 1 more detail)  ====
% (try changing these params until basis can accurately represent STA).
ktbasprs.neye = 0; % number of "identity"-like basis vectors
ktbasprs.ncos = 8; % number of raised cosine basis vectors
ktbasprs.kpeaks = [0 nkt/2+3]; % location of 1st and last basis vector bump
ktbasprs.b = 7; % determines how nonlinearly to stretch basis (higher => more linear)
[ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, nkt); % make basis
filtprs_basis = (ktbas'*ktbas)\(ktbas'*sta);  % filter represented in new basis
sta_basis = ktbas*filtprs_basis;

% Insert filter basis into fitting struct
pp0.k = sta_basis; % insert sta filter
pp0.kt = filtprs_basis; % filter coefficients (in temporal basis)
pp0.ktbas = ktbas; % temporal basis
pp0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

%% == 4. ML/MID 1:  ML estimator for LNP with CBF (cylindrical basis func) nonlinearity

nfilts = 2; % number of filters to recover

% Set parameters for cylindrical basis funcs (CBFs) and initialize fit
fstructCBF.nfuncs = 5; % number of basis functions for nonlinearity
fstructCBF.epprob = [.01, 0.99]; % cumulative probability outside outermost basis function peaks
fstructCBF.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function
fstructCBF.nlfuntype = 'cbf';

% Fit the model (iteratively adding one filter at a time)
[ppcbf,negLcbf] = fitLNP_multifilts_ML(pp0,Stim_tr,sps_tr,nfilts,fstructCBF,istacFilts);

% Compute test log-likelihood
LLcbf_tr = logli_LNP(ppcbf,Stim_tr,sps_tr); % training log-likelihood
LLcbf_tst = logli_LNP(ppcbf,Stim_tst,sps_tst); % test log-likelihood

% Compute gridded nonlinearity (for plotting purposes)
[xcbf,ycbf,nlcbf] = compNlin_2D(ppcbf.k,ppcbf.nlfun,Stim_tr); % compute gridded 2D nonlinearity

% Plot filters and nonlinearity
subplot(243);  % Plot filters
plot(tt,filts_true,'k--',tt,squeeze(ppcbf.k),'linewidth',2); 
axis tight; title('ML-cbf filters');
subplot(247); mesh(xcbf,ycbf,nlcbf); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D nonlinearity');


%% == 5. ML / MID 2:  ML estimator for LNP with RBF (radial basis func) nonlinearity

% Set parameters for cylindrical basis funcs (CBFs) and initialize fit
fstructRBF.nfuncs = 3; % number of basis functions for nonlinearity
fstructRBF.epprob = [.01, 0.99]; % cumulative probability outside outermost basis function peaks
fstructRBF.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function
fstructRBF.nlfuntype = 'rbf';

% Fit the model (iteratively adding one filter at a time)
[pprbf,negLrbf] = fitLNP_multifilts_ML(pp0,Stim_tr,sps_tr,nfilts,fstructRBF,istacFilts);

% Compute test log-likelihood
LLrbf_tr = logli_LNP(pprbf,Stim_tr,sps_tr); % training log-likelihood
LLrbf_tst = logli_LNP(pprbf,Stim_tst,sps_tst); % test log-likelihood

% Compute gridded nonlinearity (for plotting purposes)
[xrbf,yrbf,nlrbf] = compNlin_2D(pprbf.k,pprbf.nlfun,Stim_tr); % compute gridded 2D nonlinearity

% Plot filters and nonlinearity
subplot(244);  % Plot filters
plot(tt,filts_true,'k--',tt,squeeze(pprbf.k),'linewidth',2); 
axis tight; title('ML-cbf filters');
subplot(248); mesh(xrbf,yrbf,nlrbf); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D nonlinearity');



%% 6. ==== Compute training and test performance in bits/spike =====

% ==== report R2 error in reconstructing first two filters =====
ktrue = filts_true(:,1:2); % true filters 
kmse = sum((ktrue(:)-mean(ktrue(:)).^2)); % mse of these two filters around mean
ferr = @(k)(ktrue-k*(k\ktrue)); % error in optimal R2 reconstruction of ktrue
fRsq = @(k)(1-sum(sum(ferr(k).^2))/kmse); % R-squared

Rsqvals = [fRsq(istacFilts(:,1:2)),fRsq(istacFilts(:,2:3)),fRsq(squeeze(ppcbf.k)),fRsq(squeeze(pprbf.k))];

fprintf('\n=========== RESULTS =================\n');
fprintf('\nFilter R^2:\n------------\n');
fprintf('istac-12:%.2f  istac-23:%.2f  cbf:%.2f rbf:%.2f\n',Rsqvals);


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

SSinfo_tr = [f1(LListac12_tr), f1(LListac23_tr), f1(LLcbf_tr) f1(LLrbf_tr)];
SSinfo_tst = [f2(LListac12_tst),f2(LListac23_tst), f2(LLcbf_tst) f2(LLrbf_tst)];

fprintf('\nSingle-spike information (bits/spike):\n');
fprintf('------------------------------------- \n');
fprintf('Train: istac-12: %.2f  istac-23: %.2f  cbf:%.2f  rbf:%.2f\n', SSinfo_tr);
fprintf('Test:  istac-12: %.2f  istac-23: %.2f  cbf:%.2f  rbf:%.2f\n', SSinfo_tst);
