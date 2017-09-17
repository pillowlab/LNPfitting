function [pp,neglogli,filterpicked] = addfilterLNP_cbfNlin(pp,Stim,sps,kcandidates)
% [pp,neglogli,filterpicked] = addfilterLNP_cbfNlin(pp,Stim,sps,kcandidates)
%
% Picks a filter to add to an LNP model based on which gives the greatest
% increase in log-likelihood.
%
%  INPUTS: 
%             pp [1x1] -  param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%    kcandidates [1x1] - each column is a candidate filter for inclusion
%             
%          ---- Fields: ------------------------------------------
%          nfuncs = # of Gaussian basis functions
%          epprob = [low,high] quantiles for first and last ctrs
%             sig = stdev of each Gaussian function
%          nloutfun = output function (maps weighted sum of functions to positive spike rate)
%          --------------------------------------------------------
%
%  OUTPUTS:
%          ppnew [1x1] - new param struct (with estimated params)
%   negloglivals [1xP] - negative log-likelihood at ML estimate of
%                        nonlinearity, for each filter
%
% ----------------
% Algorithm:
%  (1) convolves the filters in pp.k with the stimulus
%  (2) sets the centers of the cbf Gaussian basis functions (using 1st filter only)
%  (3) fits the weights {w_i} by ML, with spike rate given by
%       lambda = g(w11*f1(x1) + w21*f2(x1) + ... wm1*fm(x1) + ...
%                  w12*f1(x2) + w22*f2(x2) + ... wm2*fm(x2) + ...
%                  ...
%                  w1n*f1(xn) + w2n*f2(x2) + ... wmn*fm(xn) )
%   where g = output function 'fstruct.nloutfun' (e.g., @exp)
%        fj = j'th basis function (out of m total basis funct)
%        xi = output of i'th filter (out of n total filters)
%
%  updated: 24 Apr, 2012 (JW Pillow)
 
% Check inputs
RefreshRate = pp.RefreshRate; % Stimulus frame rate (frames / sec)
kk = pp.k;  % filters in current model
[nkt,nkx,nfilts] = size(kk); % size and number of filters in LNP model
nkprs = nkt*nkx; % number of parameters in a single filter
kk = reshape(kk,nkprs,nfilts);
[len,nkcan] = size(kcandidates);
if len~=nkprs
    error('candidate filters don''t match size of current filters');
end

%% ==== Pre-processing of candidate filters ====================

% Project candidate filters into temporal basis ktbas:
ktbas = pp.ktbas;
kcandidates = reshape(kcandidates,nkt,nkx,nkcan);
for j = 1:nkcan
    kcandidates(:,:,j) = ktbas*(ktbas\kcandidates(:,:,j));
end
kcandidates = reshape(kcandidates,nkprs,nkcan);

% Determine which candidate filters are linearly dependent with current filters
kcandidates = normalizecols(kcandidates); % make unit vectors
kcanHat = kk*(kk\kcandidates); % best reconstruction
isLinDep = sum((kcanHat-kcandidates).^2,1)<0.001; 
if all(isLinDep)
    error('All candidate filters are linearly dependent with the current filter set');
end
if any(isLinDep)
    fprintf('Filter(s) in span of current filters: %s \n', num2str(find(isLinDep)));
end
knumskeep = find(~isLinDep);
kcandidates = kcandidates(:,~isLinDep);
nkcan = size(kcandidates,2);

% Orthogonalize candidate filters w.r.t. current filters
for j = 1:nkcan
    kcandidates(:,j) = gsorth(kcandidates(:,j),kk);
end


%% ==== Compute existing filter outputs =============================
slen = size(Stim,1); % stimulus length
iiLi = computeMask_LNP(pp.mask,slen); % compute mask (time bins to use)

% Convolve stimulus with first filter and apply mask
slen = size(Stim,1);
xproj = zeros(slen,nfilts);
for j = 1:nfilts
    xproj(:,j) = sameconv(Stim,pp.k(:,:,j));  % filter stim with filter
end
xproj = xproj(iiLi,:); % keep only those time bins within the mask
sps = sps(iiLi); % keep spikes only within mask time bins


%% ==== Compute candidate filter outputs =============================
xprojNew = zeros(slen,nkcan);
kcandidates = reshape(kcandidates,nkt,nkx,nkcan);
for j = 1:nkcan
    xprojNew(:,j) = sameconv(Stim,kcandidates(:,:,j)); % filter stim with filter
end
xprojNew = xprojNew(iiLi,:);

%% ==== Evaluate logli at ML nonlinearity for each filter ============
dtbin = 1./RefreshRate;
fprsCurrent = pp.fprs;
fprsNew = ones(pp.fstruct.nfuncs,1)*.1;
fprs0 = [fprsCurrent;fprsNew];  % initial estimate at nonlinearity params
nlfun = cell(nkcan,1); fprs=cell(nkcan,1); negL = zeros(nkcan,1);
for j = 1:nkcan
    % Fit function by maximum likelihood (Newton's method)
    [nlfun{j},fprs{j},negL(j)] = fit_nonparF_LNPmodel([xproj,xprojNew(:,j)],sps,pp.fstruct,dtbin,fprs0);
end

[neglogli,jmin] = min(negL);
filterpicked = knumskeep(jmin);

pp.nlfun = nlfun{jmin};
pp.fprs = fprs{jmin};
pp.kt(:,:,nfilts+1) = ktbas\kcandidates(:,:,jmin);
pp.k(:,:,nfilts+1) = kcandidates(:,:,jmin);  % should match ktbas*pp.kt(:,:,nfilts+1);


% =========================================
function [nlfun,fprs,neglogli] = fit_nonparF_LNPmodel(x,y,fstruct,dtbin,prs0)
% [nlfun,fwts,neglogli] = fit_nonparF_LNPmodel(x,y,fstruct,dtbin,prs0)
% 
% [nlfun,fwts,neglogli] = fit_nonparF_LNPmodel(x,y,fstruct,dtbin,prs0)
% 
% Fit nonparametric nonlinearity in LNP model in a basis of Gaussian
% functions
%
% Inputs:
%       x [Nx1] - input variable
%       y [Nx1] - output (count) variable
%   fstruct [struct] - struct with params for basis functions
%   dtbin [1x1]- time bin size (OPTIONAL; afstructumed 1 otherwise)
%   prs0 [Mx1] - initial guess at spline params (OPTIONAL)
% 
% Outputs:
%   nlfun - function handle for nonlinearity (uses 'ppval')
%   fwts  - estimated weights on basis functions
%   neglogli - negative loglikelihood at parameter estimate


if nargin<4
    dtbin=1; % size of a time bin (in secs)
end

% Compute design matrix
Xdesign = evalCBFs(x,fstruct);
nprs = size(Xdesign,2);

% Set initial params, if necessary
if nargin<5
    prs0 = randn(nprs,1)*.1;
end

% Set up loss function
lossfun = @(prs)neglogli_LNP_linearBasis(prs,Xdesign,y,fstruct.nloutfun,dtbin);

% minimize negative log-likelihood using Newton's method
[fprs,neglogli] = fminNewton(lossfun,prs0);  

% Create function handle for resulting tspline
nlfun = @(x)evalCBFnlin(x,fstruct,fprs);

% -------------------
% for debupping only
% -------------------
% HessCheck(lossfun,prs0); % Check that analytic grad and Hessian are correct
