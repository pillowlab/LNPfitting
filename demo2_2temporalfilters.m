% demo2_2temporalfilters.m
%
% Tutorial script illustrating maximum likelihood / maximally informative
% dimensions (MID) estimation for LNP model with TWO temporal filters

% initialize paths
initpaths;

datasetnum = 2;  % select: 1 (white noise) or 2 (correlated)
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
nFilts = 3; % number of filters to compute
[k_istac,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nFilts); % find iSTAC filters

% Fit iSTAC model nonlinearity using filters 1 and 2
pp_istac12 = fitNlin_expquad_ML(Stim_tr,sps_tr,k_istac(:,1:2),RefreshRate); % LNP model struct

% Visualize 2D nonlinearity
[xgrd,ygrd,nlvals] = compNlin_2D(pp_istac12.k,pp_istac12.nlfun,Stim_tr); % compute gridded version
subplot(221); mesh(xgrd,ygrd,nlvals); 
zlm = [0 max(nlvals(:))*1.01];
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 1 & 2)')

% Fit iSTAC model nonlinearity using filters 2 and 3
pp_istac23 = fitNlin_expquad_ML(Stim_tr,sps_tr,k_istac(:,2:3),RefreshRate); % LNP model struct

% Visualize 2D nonlinearity
[xgrd,ygrd,nlvals] = compNlin_2D(pp_istac23.k,pp_istac23.nlfun,Stim_tr);
subplot(222); mesh(xgrd,ygrd,nlvals); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 2 axis');ylabel('filter 3 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 2 & 3)')


%% == 3. Set up temporal basis for stimulus filters (for ML / MID estimators)

pp0 = makeFittingStruct_LNP(k_istac(:,1),RefreshRate); % initialize param struct

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

%% == 4. MID1:  ML estimator for LNP with CBF (cylindrical basis func) nonlinearity

% Write single func to do this (up to N filters)


%% == 5. MID2:  ML estimator for LNP with RBF (radial basis func) nonlinearity

% Write single func to do this (up to 3 filters)
