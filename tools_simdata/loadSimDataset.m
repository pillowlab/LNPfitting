function [Stim_tr,sps_tr,Stim_tst,sps_tst,RefreshRate,filts_true] = loadSimDataset(datasetnum,trainfrac)
% [Stim_tr,sps_tr,Stim_tst,sps_tst,RefreshRate,filts_true] = loadSimDataset(datasetnum)
%
% Loads (and/or creates) simulated datasets for LNPfitting demo code
%
% INPUT
%  datasetnum - 1, temporal white noise
%               2, temporal correlated noise
%
% OUTPUT
%    Stim_tr [ntrain x nx] - training stimulus
%     sps_tr [ntrain x 1]  - training spike train (counts / bin)
%   Stim_tst [ntest x nx]  - test stimulus
%    sps_tst [ntest x 1]   - test spike train (counts / bin)
%    RefreshRate [1 x 1]   - assumed frame rate of stimulus (frames / sec)
%             filts_true   - true filters of simulated model

% pick dataset to load (or create if necessary)
switch datasetnum
    case 1
        datasetname = 'simdatadir/simdata1.mat';  % name of dataset
        if ~exist(datasetname,'file') % Create simulated dataset if necessary
            fprintf('Creating simulated dataset: ''%s''\n', datasetname);
            mkSimData1_1Dwhitenoisestim;
        end
    case 2
        datasetname = 'simdatadir/simdata2.mat';  % name of dataset
        if ~exist(datasetname,'file') % Create simulated dataset if necessary
            fprintf('Creating simulated dataset: ''%s''\n', datasetname);
            mkSimData2_1Dcorrstim;
        end
end

% Load data
load(datasetname); % load dataset (loads struct 'simdata')
slen = size(simdata.Stim,1); % number of time bins in full stimulus / spike train
slen_tr = round(trainfrac*slen); % number of time bins in training dataset
    
% Set training data
Stim_tr = simdata.Stim(1:slen_tr,:);
sps_tr = simdata.spikes(1:slen_tr,:);

% Set test data
Stim_tst = simdata.Stim(slen_tr+1:end,:);
sps_tst = simdata.spikes(slen_tr+1:end,:);

% Load other variables
RefreshRate = simdata.RefreshRate; % stimulus refresh rate (in Hz).
filts_true = simdata.filts_true; % true LNP model filters

% Report length of training data and spike counts / rates
nsp_tr = sum(sps_tr); % Determine how many spikes in training set
fprintf('\n------------\nLoaded %s\n',datasetname);
fprintf('Total length: %d bins (training data: %d bins)\n', slen, slen_tr);
fprintf('Number of spikes in training data: %d (%.2f sp/sec)\n', nsp_tr, nsp_tr/slen_tr*RefreshRate);

