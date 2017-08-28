% tutorial1_MID_tfilter.m
%
% Tutorial script for illustrating MID estimator code for neuron with
% purely temporal stimulus

% initialize paths
initpaths;

% 0. Load dataset (from simulated LNP neuron)
if ~exist('simdata/simdata1.mat','file')
    fprintf('Creating simulated dataset ''simdata1.mat''\n');
    mkSimDataset1_tfilterLNP;
end

load simdata/simdata1.mat;  % load dataset
slen = size(simdata1.Stim,1);

nkt = 30; % number of time bins in filter
% (Normally would want to experiment to set this (eg by inspecting the STA) 


%% 1. Divide into training and test datasets

trainfrac = .8; % fraction of data to set aside as "test data"
slen_tr = round(trainfrac*slen); % length of training dataset
slen_test = slen-slen;  % length of test dataset
    
% Set training data
Stim_tr = simdata1.Stim(1:slen_tr,:);
sps_tr = simdata1.spikes(1:slen_tr,:);

% Set test data
Stim_test = simdata1.Stim(slen_tr+1:end,:);
sps_test = simdata1.spikes(slen_tr+1:end,:);

nsp = sum(sps_tr);
fprintf('Number of spikes in training dataset: %d\n', nsp);
    

%% 2. Initialize filters using iSTAC estimator

% Compute STA and STC
fprintf('Computing STA and STC...\n');
[sta,stc,rawmu,rawcov] = simpleSTC(Stim_tr,sps_tr,nkt);

% Compute first iSTAC filter 
nfilters = 3;  % 
condthresh = 0.05; % threshold on condition number (for pruning low-variance stimulus axes)
fprintf('Computing top iSTAC filter...\n');
filts_init = compiSTAC(sta,stc,rawmu,rawcov,nfilters,condthresh);

% Initialize struct for fitting 1-filter model
mask = [];  % use all training data
gg0 = makeFittingStruct_LNP(filts_init(:,1),RefreshRate,mask); % create LNP fitting structure


%% 3. Set up temporal basis for representing MID filters

% Set parameters for temporal basis and inspect accuracy of reconstruction
ktbasprs.neye = 0; % number of "identity"-like basis vectors
ktbasprs.ncos = 10; % number of raised cosine basis vectors
ktbasprs.kpeaks = [1 nkt/2+3]; % location of 1st and last basis vector bump
ktbasprs.b = 1.5; % determines how nonlinearly to stretch basis (higher => more linear)
[ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, nkt); % make basis
filtprs_basis = (ktbas\filts_init(:,1));  % filter represented in new basis
filt_basis = ktbas*filtprs_basis;

% Insert filter into new fitting struct
gg0.k = filt_basis; % filter
gg0.kt = filtprs_basis; % filter coefficients (in temporal basis)
gg0.ktbas = ktbas; % temporal basis
gg0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

% Plot iSTAC vs. best reconstruction in temporal basis
ttk = -nkt+1:0;
subplot(211); % ----
plot(ttk,ktbasis); xlabel('time bin'); title('temporal basis'); axis tight;
subplot(212); % ----
plot(ttk,filts_init(:,1),'b',ttk,filt_basis(:,1),'r--'); axis tight; title('iSTAC filter vs. basis fit');


%% 4.  Now estimate filters using MID estimator with RBF nonlinearity

fstruct.nfuncs = 3; % number of basis functions for nonlinearity
fstruct.epprob = [0, 1]; % cumulative probability outside outermost basis function peaks
fstruct.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function

fprintf('Computing MID estimator, 1-filter model\n');
%[gg0,negL0_tr] = fitNlin_CBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate nonparametric nonlinearity
[gg0r,negL0r_tr] = fitNlin_RBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate RBF nonlinearity (should be identical)
    

