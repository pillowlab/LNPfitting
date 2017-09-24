function [pp,negL,pp_prev,negL_prev] = fitLNP_multifiltsCBF_ML(pp,Stim,sps,nfilts,fstruct,initFilts,optimArgs)
% [pp,negL,pp_prev,negL_prev] = fitLNP_multifiltsCBF_ML(pp,Stim,sps,nfilts,fstruct,initFilts,optimArgs)
%
% Maximum likelihood / MID fitting of LNP model with multiple filters and
% either CBF parametrized nonlinearity
%
%  INPUTS:
%             pp [1x1] - param struct for LNP model
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%         nfilts [1x1] - max number of filters to fit
%  fstruct [structure] - nonlinearity structure with fields
%                        .nlfuntype ('rbf' or 'cbf')
%                        .nfuncs (# of basis functions)
%                        .epprob (quantiles for 1st and last basis func, e.g. [0.01, 0.99])
%                        .nloutfun (final output nonlinearity)
%   initFilts [NM x P] - candidate filters for initializing fits (optional), where P >= nfilts
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  OUTPUS:
%         ppnew [1x1] - new param struct (with estimated params)
%          negL [1x1] - negative log-likelihood at ML estimate
%       pp_prev [L-1x1] - cell array of model structs with 1 to L-1 filters
%                         (obtained en route to fitting full model) 
%     negL_prev [L-1x1] - array of negative log-likelihood of smaller models
%                         (with 1 to L-1 filters) on the training data  


% Check nonlinearity type
if ~strcmpi(fstruct.nlfuntype,'cbf')
        error('wrong type of nonlinearity: fstruct.nlfuntype should be ''cbf''');
else
    pp.nlfuntype = fstruct.nlfuntype;
end

% Extract stimulus temporal basis and filter size
ktbas = pp.ktbas;
nkt = size(ktbas,1); % number of time bins in filter
nkx = size(Stim,2);  % number of spatial pixels in filter

% Set optimization arguments if necessary
if nargin<7
    optimArgs = {'display','on'};
end

%% Compute initial filters using iSTAC if none provided
if nargin<6 || isempty(initFilts)
    fprintf('\nRunning iSTAC to get initial filter estimates...\n');
    [sta,stc,rawmu,rawcov] = simpleSTC(Stim,sps,nkt);  % compute raw and spike-trippered moments
    initFilts = compiSTAC(sta(:),stc,rawmu,rawcov,nfilts); % compute iSTAC filters
elseif size(initFilts,2)<nfilts
    warning('fewer filters in ''initFilts'' than # filters requested; aumenting with iSTAC filters');
    % now run iSTAC to supplement filters
    fprintf('\nRunning iSTAC to supplement initial filter estimates...\n');    
    [sta,stc,rawmu,rawcov] = simpleSTC(Stim,sps,nkt);  % compute raw and spike-trippered moments
    initFilts0 = compiSTAC(sta(:),stc,rawmu,rawcov,nfilts); % compute iSTAC filters
    initFilts = [initFilts, initFilts0];
end

%% Insert first filter into param struct as represented in basis
filt0full = reshape(initFilts(:,1),nkt,nkx); % full initial filter
pp.kt = (ktbas\filt0full); % temporal basis coeffs for filter
pp.k = ktbas*pp.kt; % full filter represented in temporal basis


%% Fit 1-filter model
fprintf('\n================\nfitLNP_multifiltsCBF_ML: Fitting filter 1 (of %d)\n================\n',nfilts);
[pp0,~] = fitNlin_CBFs(pp,Stim,sps,fstruct); % Initialize estimate of nonlinearity
[pp,negL] = fitLNP_multifilts_cbfNlin(pp0,Stim,sps,optimArgs);  % Do joint fitting of nonlinearity and filter


%% Fit 2 and above filter models.

% Initialize variables for smaller models
pp_prev = cell(nfilts-1,1); % models with fewer filter
negL_prev = zeros(nfilts-1,1); % training negative log-likelihood

for jj = 2:nfilts

    % Store previous param struct and neg logli value
    pp_prev{jj-1} = pp;
    negL_prev(jj-1) = negL;

    % Add filter to model
    fprintf('\n================\nfitLNP_multifiltsCBF_ML: Fitting filter %d (of %d)\n================\n',jj,nfilts);

    [pp0,~,filterPicked] = addfilterLNP_cbfNlin(pp,Stim,sps,initFilts); % pick a filter to add and initialize nonlinearity

    % Now optimize negative log-likelihood for filter and nonlinearity
    fprintf('\nInitializing with istac filter #%d\n',filterPicked);
    [pp,negL] = fitLNP_multifilts_cbfNlin(pp0,Stim,sps,optimArgs); % jointly fit filter and nonlinearity
    
end

% store final object in same cell array (for convenience)
pp_prev{nfilts} = pp;
negL_prev(nfilts) = negL;
