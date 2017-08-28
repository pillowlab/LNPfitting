% tutorial1_MID_tfilter.m
%
% Tutorial script for illustrating MID estimator code for neuron with
% purely temporal stimulus

% 0. Load dataset (from simulated LNP neuron)
if ~exist('simdata/simdata1.mat','file')
    fprintf('Creating simulated dataset ''simdata1.mat''\n');
    mkSimDataset1_tfilterLNP;
end

slen = size(simdata1.Stim,1);

%% 1. Divide into training and test datasets

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
    

%% 2. Initialize filters using iSTAC estimator

% Compute STA and STC
fprintf('Computing STA and STC...\n');
[sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes,nkt);

% Compute first iSTAC filter 
nfilters = 3;  % 
condthresh = 0.05; % threshold on condition number (for pruning low-variance stimulus axes)
fprintf('Computing top iSTAC filter...\n');
filts_init = compiSTAC(sta,stc,rawmu,rawcov,nfilters,condthresh);

% Initialize struct for fitting 1-filter model
mask = [];  % use all training data
gg0 = makeFittingStruct_LNP(filts_init(:,1),RefreshRate,mask); % create LNP fitting structure


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
ttk = 
subplot(211); % ----
plot(ttk,ktbasis); xlabel('time bin'); title('temporal basis'); axis tight;
subplot(212); % ----
plot(ttk,filt0,'b',ttk,filt_basis,'r--'); axis tight; title('iSTAC filter vs. basis fit');


%% 5.  Now estimate filters using MID estimator with RBF nonlinearity

fstruct.nfuncs = 3; % number of basis functions for nonlinearity
fstruct.epprob = [0, 1]; % cumulative probability outside outermost basis function peaks
fstruct.nloutfun = @logexp1;  % log(1+exp(x)) % nonlinear stretching function

fprintf('Computing MID estimator, 1-filter model\n');
%[gg0,negL0_tr] = fitNlin_CBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate nonparametric nonlinearity
[gg0r,negL0r_tr] = fitNlin_RBFs(gg0,Stim_tr,sps_tr,fstruct); % estimate RBF nonlinearity (should be identical)
    

