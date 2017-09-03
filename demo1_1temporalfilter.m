% demo1_1temporalfilter.m
%
% Tutorial script for illustrating MID estimator code for a single-filter model
% purely temporal stimulus

% initialize paths
initpaths;

% Create simulated dataset if necessary
datasetname = 'simdata/simdata1.mat';  % name of dataset
if ~exist(datasetname,'file')
    fprintf('Creating simulated dataset: ''%s''\n', datasetname);
    mkSimDataset1_tfilterLNP;
end

%% 1. Load data and divide into training and tst datasets

trainfrac = .8; % fraction of data to set aside as "tst data"

% Load data
load(datasetname); % load dataset 
RefreshRate = simdata1.RefreshRate; % stimulus refresh rate (in Hz).
slen = size(simdata1.Stim,1); % number of time bins in stimulus
slen_tr = round(trainfrac*slen); % length of training dataset
slen_tst = slen-slen_tr;  % length of tst dataset
    
% Set training data
Stim_tr = simdata1.Stim(1:slen_tr,:);
sps_tr = simdata1.spikes(1:slen_tr,:);
% Set test data
Stim_tst = simdata1.Stim(slen_tr+1:end,:);
sps_tst = simdata1.spikes(slen_tr+1:end,:);

nsp_tr = sum(sps_tr); % Determine how many spikes in training set
fprintf('\n------------\nLoaded %s\n',datasetname);
fprintf('Total length: %d bins (training data: %d bins)\n', slen, slen_tr);
fprintf('Number of spikes in training data: %d (%.2f sp/sec)\n', nsp_tr, nsp_tr/slen_tr*RefreshRate);


%% == 2. Compute STA and estimate (piecewise constant) nonlinearity using histograms ====

nkt = 30; % number of time bins to use for filter 
% This is important: normally would want to vary this to see how many bins
% we need to capture structure in the STA

% Compute STA
sta = simpleSTC(Stim_tr,sps_tr,nkt);  % compute STA
sta = sta./norm(sta);  % normalize sta to be a unit vector

% Estimate piecewise-constant nonlinearity ("reconstruction" method using histograms)
nhistbins = 15; % # histogram bins to use
fnlhist = fitNlin_hist1D(Stim_tr, sps_tr, sta, RefreshRate, nhistbins); % estimate 1D nonlinearity 

%% == 3.  Fit LNP model with exponential nonlinearity  ================

slen = size(Stim_tr,1);  % total number of time bins in stimulus 
mask = [];  % time range to use for fitting (set to [] if not needed).

% Set up fitting structure and compute initial logli
pp0 = makeFittingStruct_LNP(sta,RefreshRate,mask); % param struct
negL0 = -logli_LNP(pp0,Stim_tr,sps_tr);  % negative log-likelihood
fprintf('\nFitting LNP model w/ exp nonlinearity\n');

% Do ML estimation of model params (with temporal basis defined in gg0)
opts = {'display', 'off', 'maxiter', 100};
[pp_exp,negL1,C1] = fitLNP_1filt_ML(pp0,Stim_tr,sps_tr,opts); % find MLE by gradient ascent
eb1 = sqrt(diag(C1(1:nkt,1:nkt))); % 1SD error bars


%% == 4.  Run MID: estimate filter and non-parametric nonlinearity with rbf basis ===

% Set parameters for radial basis functions (RBFs), for parametrizing nonlinearity
fstruct.nfuncs = 7; % number of RBFs
fstruct.epprob = [.01, .99]; % cumulative probability outside outermost basis function peaks (endpoints)
fstruct.nloutfun = @logexp1;  % log(1+exp(x))  % nonlinear output function
fprintf('\nFitting LNP model w/ rbf nonlinearity\n');
 
% Do fitting
[pp_rbf,negLnonpar0] = fitNlin_CBFs(pp_exp,Stim_tr,sps_tr,fstruct);  % initialize nonlinearity while holding filter fixed
opts = {'display', 'off'}; % optimization parameters
[pp_rbf,negLnonpar] = fitLNP_multifilts_cbfNlin(pp_rbf,Stim_tr,sps_tr,opts); % jointly fit filter and nonlinearity


%% 5. ====== Make plots & report performance ============

% -- plot filters (rescaled as unit vectors) ---
tt = (-nkt+1:0)/RefreshRate;  % time points in Stim_tr filter

% true 1st filter (from simulation)
uvec = @(x)(x./norm(x)); % anonymous function to create unit vector
trueK = uvec(simdata1.filts_true(:,1));

% -- Plot filter and filter estimates ---------
subplot(221); 
plot(tt,trueK,'k',tt,sta,tt,uvec(pp_exp.k),tt,uvec(pp_rbf.k), 'linewidth',2);
legend('true','sta','ML-exptl','ML-rbf','location', 'northwest');
xlabel('time before spike (ms)'); ylabel('weight');
title('filters (rescaled as unit vectors)');  axis tight;

% ---- Compute nonlinearities for plotting ---------
xnl = (min(Stim_tr)-.25):.1:(max(Stim_tr)+.25); % x points for evaluating nonlinearity
ynl_hist = fnlhist(xnl); % histogram-based (piecewise constant) nonlinearity 
ynl_exp = exp(xnl*norm(pp_exp.k)+pp_exp.dc);  % exponential nonlinearity
ynl_rbf = pp_rbf.nlfun(xnl*norm(pp_rbf.k));   % rbf nonlinearity

% ---- Plot nonlinearities --------------------------
subplot(222); 
plot(xnl, ynl_hist,xnl,ynl_exp,xnl,ynl_rbf, 'linewidth',2);
axis tight; set(gca,'ylim',[0 200]);
ylabel('rate (sps/s)'); xlabel('filter output');
legend('hist','ML-exptl','ML-rbf','location', 'northwest');
title('estimated nonlinearities');

% ==== report filter estimation error =========
ferr = @(k)(1-(norm(uvec(k)-trueK)));  % normalized error
fprintf('\nFilter R^2:\n------------\n');
fprintf('sta:%.2f  exp:%.2f  rbf:%.2f\n', [ferr(sta), ferr(pp_exp.k), ferr(pp_rbf.k)]);

% ==== Compute training and test performance in bits/spike =====

% Compute the log-likelihood under constant rate (homogeneous Poisson) model
nsp_tst = sum(sps_tst);  % number of spikes in tst set
muspike_tr = nsp_tr/slen_tr;       % mean number of spikes / bin, training set
muspike_tst = nsp_tst/slen_tst; % mean number of spikes / bin, tst set
LL0_tr =   nsp_tr*log(muspike_tr) - slen_tr*muspike_tr; % log-likelihood, training data
LL0_tst = nsp_tst*log(muspike_tst) - slen_tst*muspike_tst; % log-likelihood test data

% 1. Compute logli for lnp with histogram nonlinearity
pp_sta = pp0; % make struct for the sta+histogram-nonlinearity model
pp_sta.k = sta; % insert STA as filter
pp_sta.dc = 0; % remove DC component (if necessary)
pp_sta.kt = []; pp_sta.ktbas = []; % remove basis stuff (just to make sure it isn't used accidentally)
pp_sta.nlfun =  fnlhist;
LLsta_tr = logli_LNP(pp_sta,Stim_tr,sps_tr); % training log-likelihood
[LLsta_tst,rrsta_tst] = logli_LNP(pp_sta,Stim_tst,sps_tst); % test log-likelihood

% 2. Compute logli for lnp with exponential nonlinearity
%ratepred_pGLM = exp(pGLMconst + Xdsgn*pGLMfilt); % rate under exp nonlinearity
LLexp_tr = logli_LNP(pp_exp,Stim_tr,sps_tr); % train
[LLexp_tst,rrexp_tst] = logli_LNP(pp_exp,Stim_tst,sps_tst); % test

% 3. Compute logli for lnp with rbf nonlinearity
LLrbf_tr = logli_LNP(pp_rbf,Stim_tr,sps_tr); % train
[LLrbf_tst,rrrbf_tst] = logli_LNP(pp_rbf,Stim_tst,sps_tst); % test

% Single-spike information:
% ------------------------
% The difference of the loglikelihood and homogeneous-Poisson
% loglikelihood, normalized by the number of spikes, gives us an intuitive
% way to compare log-likelihoods in units of bits / spike.  This is a
% quantity known as the ("empirical" or "sample") single-spike information.
% [See Brenner et al, Neural Comp 2000; Williamson et al, PLoS CB 2015].

f1 = @(x)((x-LL0_tr)/nsp_tr/log(2)); % compute training single-spike info 
f2 = @(x)((x-LL0_tst)/nsp_tst/log(2)); % compute test single-spike info
% (if we don't divide by log 2 we get it in nats)

SSinfo_tr = [f1(LLsta_tr), f1(LLexp_tr), f1(LLrbf_tr)];
SSinfo_tst = [f2(LLsta_tst), f2(LLexp_tst), f2(LLrbf_tst)];

fprintf('\nSingle-spike information (bits/spike):\n');
fprintf('------------------------------------- \n');
fprintf('Train: sta-hist:%.2f  exp:%.2f  rbf:%.2f\n', SSinfo_tr);
fprintf('Test:  sta-hist:%.2f  exp:%.2f  rbf:%.2f\n', SSinfo_tst);

% ==== Last: plot the rate predictions for the two models =========
subplot(212); 
iiplot = 1:200; % time bins to plot
stem(iiplot,sps_tst(iiplot), 'k'); hold on;
plot(iiplot,rrsta_tst(iiplot)/RefreshRate, ...
    iiplot,rrexp_tst(iiplot)/RefreshRate, ...
    iiplot,rrrbf_tst(iiplot)/RefreshRate,'linewidth',2); 
 hold off; title('rate predictions on test data');
ylabel('spikes / bin'); xlabel('time (bins)');
set(gca,'xlim', iiplot([1 end]));
legend('spike count', 'sta-hist', 'ML-exp', 'ML-rbf');
