function gg = makeFittingStruct_LNP(k0,RefreshRate,mask)
% makeFITTINGstruct_LNP - make a fitting structure for an LNP model
%
% gg = makeFittingStruct_LNP(k0,RefreshRate,mask)
%
% Inputs:  
%         k0 [nkt x nkx] - initial guess at kernel (can pass all zeros)
%           OR
%        gg0 [struct] - structure with params gg0.k and gg0.dc }
% RefreshRate [1 x 1] - stimulus refresh rate (frames per second)
%        mask [M x 2] - each row specify time bins to use for fitting (OPTIONAL)
%
% Ouptuts:
%        gg [struct] - default structure for lnp fitting
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
    gg = k0;
else
    gg.k = k0;  %  stim filter k
    gg.dc = 0;
    gg.kt = [];
    gg.ktbas = [];
    gg.ktbasprs = [];
end

% set nonlinearity
if ~isfield(gg, 'nlfun') || isempty(gg.nlfun)
    gg.nlfun = @expfun;
end

% ==================================================================
% Set up default temporal basis for stim kernel
% ==================================================================
if ~isfield(gg, 'ktbas') || isempty(gg.ktbas)
    nkt = size(gg.k,1);  % number of temporal elements in the k
    ktbasprs.neye = min(5,floor(nkt/2)); % # "identity" basis vectors near time of spike;
    ktbasprs.ncos = min(5,floor(nkt/2)); % # raised-cosine vectors to use
    ktbasprs.kpeaks = [0 ((nkt-ktbasprs.neye)/2)]; % Position of 1st and last bump
    ktbasprs.b = 1; % Offset for nonlinear scaling (larger -> more linear)
    ktbas = makeBasis_StimKernel(ktbasprs,nkt);
    gg.ktbas = ktbas;
    gg.ktbasprs = ktbasprs;
end

% ==================================================================
% set up initial K params (stim filter, linearly parametrized)
% ==================================================================
gg.kt = gg.ktbas\gg.k;  % least-sqaures fit of k in temporal basis 
gg.k = gg.ktbas*gg.kt;

gg.RefreshRate = RefreshRate;
gg.mask = mask;
gg.model = 'LNP';  % 
gg.ktype = 'linear';  % type of filter parameterization (can be 'linear' or 'bilinear')
