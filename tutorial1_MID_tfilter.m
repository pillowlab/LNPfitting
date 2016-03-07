% tutorial1_MID_tfilter.m
%
% Tutorial script for illustrating MID estimator code for neuron with
% purely temporal stimulus

addpath code_LNPfitting
addpath nlfuns;
addpath code_iSTAC

% 0. Set up Linear-Nonlinear-Poisson (LNP) model neuron 
nkt = 32;           % number of time bins in filter
tk = (-nkt+1:0)'; % vector of time indices (in units of stim frames)
dtSim = .01; % bin size for representing time (here => 100 Hz frame rate)
ttk = tk*dtSim;

% Make some fake filters 
filt1 = exp(-((tk+nkt/4)/(nkt/10)).^2)-.25*exp(-((tk+nkt/2)/(nkt/4)).^2); % difference of Gaussians
filt1 = filt1./norm(filt1);  %normalize

filt2 = [diff(filt1); 0];  % 2nd filter
filt2 = filt2./norm(filt2); % normalize

filt3 = [diff(filt2); 0];  % 3rd filter
filt3 = filt3./norm(filt3); % normalize

% Plot these filters
plot(ttk, [filt1 filt2 filt3]);  
title('filters for simulation');
xlabel('time before spike (s)'); ylabel('filter coeff');
axis tight;

%% 1.  Simulate data from LNP neuron

% Create stimulus 
slen = 20000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);
Stim = conv2(Stim,normpdf(-3:3,0,1)','same'); % smooth stimulus
RefreshRate = 100; % refresh rate

% Convolve stimulus with filters
f1 = sameconv(Stim,filt1);
f2 = sameconv(Stim,filt2);
f3 = sameconv(Stim,filt3);

% Compute output of nonlinearity
softrect = @(x)(log(1+exp(x))); % soft-rectification function
fnlin = @(x1,x2,x3)(softrect(120./(1+exp(x1-1))+10*x2.^2+5*(x3).^2-80));
lam = fnlin(f1,f2,f3);

%  Simulate spike train
dtSim = 100; % in Hz
spikes = poissrnd(lam/RefreshRate); % generate spikes


%% 2. Divide into training and test datasets

trainfrac = .8; % fraction of data to set aside as "test data"
slen_tr = round(trainfrac*slen); % length of training dataset
slen_test = slen-slen;  % length of test dataset
    
% Set training data
Stim_tr = Stim(1:slen_tr,:);
sps_tr = spikes(1:slen_tr,:);

% Set test data
Stim_test = Stim(slen_tr+1:end,:);
sps_test = spikes(slen_tr+1:end,:);

nsp = sum(sps_tr);
fprintf('Number of spikes in training dataset: %d\n', nsp);
    

%% 3. Compute first iSTAC filter

% Compute STA and STC
fprintf('Computing STA and STC...\n');
[sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes,nkt);

% Compute first iSTAC filter 
nfilters = 1;  % 
condthresh = 0.05; % threshold on condition number (for pruning low-variance stimulus axes)
fprintf('Computing first iSTAC filter...\n');
filt0 = compiSTAC(sta,stc,rawmu,rawcov,nfilters,condthresh);

% Initialize struct for MID fitting of LNP model
mask = [];  % use all training data
gg0 = makeFittingStruct_LNP(filt0,RefreshRate,mask); % create LNP fitting structure


%% 4. Set up temporal basis for representing MID filters

% Set parameters for temporal basis and inspect accuracy of reconstruction
ktbasprs.neye = 0; % number of "identity"-like basis vectors
ktbasprs.ncos = 10; % number of raised cosine basis vectors
ktbasprs.kpeaks = [0 nkt/2+3]; % location of 1st and last basis vector bump
ktbasprs.b = 2.5; % determines how nonlinearly to stretch basis (higher => more linear)
[ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, nkt); % make basis
filtprs_basis = (ktbas\filt0);  % filter represented in new basis
filt_basis = ktbas*filtprs_basis;

% Insert filter into new fitting struct
gg0.k = filt_basis; % filter
gg0.kt = filtprs_basis; % filter coefficients (in temporal basis)
gg0.ktbas = ktbas; % temporal basis
gg0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

% Plot iSTAC vs. best reconstruction in temporal basis
subplot(211); % ----
plot(ttk,ktbasis); xlabel('time bin'); title('temporal basis'); axis tight;
subplot(212); % ----
plot(ttk,filt0,'b',ttk,filt_basis,'r--'); axis tight; title('iSTAC filter vs. basis fit');


%% 5.  Now estimate filters using MID estimator with RBF nonlinearity

fstruct.nfuncs = 3; % number of basis functions for nonlinearity
fstruct.epprob = [0, 1]; % cumulative probability outside outermost basis function peaks
fstruct.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function
%[gg0,negL0_tr] = fitNlin_CBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate nonparametric nonlinearity
[gg0r,negL0r_tr] = fitNlin_RBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate RBF nonlinearity (should be identical)
    

