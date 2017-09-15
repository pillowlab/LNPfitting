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
nFilts = 3; % number of filters to compute
[istacFilts,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nFilts); % find iSTAC filters

% Fit iSTAC model nonlinearity using filters 1 and 2
filts1 = istacFilts(:,1:2);
pp_istac1 = fitNlin_expquad_ML(Stim_tr,sps_tr,filts1,RefreshRate); % LNP model struct

% Visualize 2D nonlinearity
[xgrd,ygrd,nlvals] = compNlin_2D(filts1,pp_istac1.nlfun,Stim_tr);
subplot(221); mesh(xgrd,ygrd,nlvals); 
zlm = [0 max(nlvals(:))*1.01];
axis tight; set(gca,'zlim',zlm);
xlabel('filter 1 axis');ylabel('filter 2 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 1 & 2)')

% Fit iSTAC model nonlinearity using filters 2 and 3
filts2 = istacFilts(:,2:3);
pp_istac2 = fitNlin_expquad_ML(Stim_tr,sps_tr,filts2,RefreshRate); % LNP model struct

% Visualize 2D nonlinearity
[xgrd,ygrd,nlvals] = compNlin_2D(filts2,pp_istac2.nlfun,Stim_tr);
subplot(222); mesh(xgrd,ygrd,nlvals); 
axis tight; set(gca,'zlim',zlm);
xlabel('filter 2 axis');ylabel('filter 3 axis'); 
zlabel('spike rate (sps/sec)'); title('2D iSTAC nonlinearity (filters 2 & 3)')
