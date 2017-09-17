function pp = makeFittingStruct_LNP(k0,RefreshRate,mask)
% makeFITTINGstruct_LNP - make a fitting structure for an LNP model
%
% pp = makeFittingStruct_LNP(k0,RefreshRate,mask)
%
% Inputs:  
%         k0 [nkt x nkx] - initial guess at kernel (can pass all zeros)
%           OR
%        pp0 [struct] - structure with params pp0.k and pp0.dc }
% RefreshRate [1 x 1] - stimulus refresh rate (frames per second)
%        mask [M x 2] - each row specify time bins to use for fitting (OPTIONAL)
%
% Ouptuts:
%        pp [struct] - default structure for lnp fitting
%
% updated: Jan 16, 2014 (JW Pillow)

% set mask to empty if not passed in
if (nargin<2)  
    RefreshRate = 100;
    fprintf('Assumed refresh rate of %d\n',RefreshRate);
end
if (nargin<3)
    mask = [];
end
% ==================================================================
% Set up structure, if necessary
% ==================================================================
if isstruct(k0)
    pp = k0;
else
    pp.k = k0;  %  stim filter k
    pp.dc = 0;
    pp.kt = [];
    pp.ktbas = [];
    pp.ktbasprs = [];
end

% set nonlinearity
if ~isfield(pp, 'nlfun') || isempty(pp.nlfun)
    pp.nlfun = @expfun;
end

% find temporal and spatial dimensions of filter
[nt,nx] = size(k0); 

% ==================================================================
% Set up default temporal basis for stim kernel
% ==================================================================
if (nt == 1) % no temporal basis if only 1 time bin
    pp.ktbas = 1;
    pp.ktbasprs = 'none';
elseif ~isfield(pp, 'ktbas') || isempty(pp.ktbas)
    nkt = size(pp.k,1);  % number of temporal elements in the k
    ktbasprs.neye = min(5,floor(nkt/2)); % # "identity" basis vectors near time of spike;
    ktbasprs.ncos = min(5,floor(nkt/2)); % # raised-cosine vectors to use
    ktbasprs.kpeaks = [0 ((nkt-ktbasprs.neye)/2)]; % Position of 1st and last bump
    ktbasprs.b = 1; % Offset for nonlinear scaling (larger -> more linear)
    ktbas = makeBasis_StimKernel(ktbasprs,nkt);
    pp.ktbas = ktbas;
    pp.ktbasprs = ktbasprs;
end

% ==================================================================
% set up initial K params (stim filter, linearly parametrized)
% ==================================================================
pp.kt = pp.ktbas\pp.k;  % least-sqaures fit of k in temporal basis 
pp.k = pp.ktbas*pp.kt;

pp.RefreshRate = RefreshRate;
pp.mask = mask;
pp.model = 'LNP';  % 
pp.ktype = 'linear';  % type of filter parameterization (can be 'linear' or 'bilinear')
