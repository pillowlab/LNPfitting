function [pp,neglogli] = fitNlin_RBFs2(pp,Stim,sps)
% [pp,neglogli] = fitNlin_RBFs(pp,Stim,sps,fstruct)
%
% Fit nonlinearity in an LNP model with nonlinearity parameterized by 2D or
% 3D radial basis functions.
%
%  Inputs: 
%             pp [1x1] -  param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%        fstruct [1x1] -  structure with params for basis functions
%             
%          ---- Fields: ------------------------------------------
%          nfuncs = # of Gaussian basis functions
%          epprob = [low,high] quantiles for first and last ctrs
%             sig = stdev of each Gaussian function
%          nloutfun = output nonlinearity (maps weighted sum of functions to positive spike rate)
%          --------------------------------------------------------
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  Outputs:
%     ppnew [1x1] - new param struct (with estimated params)
%  neglogli [1x1] - negative log-likelihood at ML estimate
%
% updated: Jan 31, 2014 (JW Pillow)

RefreshRate = pp.RefreshRate; % stimulus frame rate
dtbin = 1./RefreshRate; % length of a single time bin
nfilts = size(pp.k,3); % number of filters in LNP model
if nfilts > 4
    error('RBFs not implemented for > 4 dimensions');
end    

% Compute mask times (times when likelihood will be computed)
slen = size(Stim,1); % stimulus length
iiLi = computeMask_LNP(pp.mask,slen); % compute mask (time bins to use)

% Convolve stimulus with filter and apply mask
slen = size(Stim,1);
Istm = zeros(slen,nfilts);
for j = 1:nfilts
    Istm(:,j) = sameconv(Stim,pp.k(:,:,j));  % filter stim with filter
end
Istm = Istm(iiLi,:); % keep only those time bins within the mask
sps = sps(iiLi); % keep spikes only within mask time bins

% Set RBF centers and sigma along each input dimension
fstruct = pp.fstruct;
while size(fstruct.ctrs,1) < nfilts
    fstruct.ctrs(end+1,:) = fstruct.ctrs(end,:);
    fstruct.sig(end+1,1) = fstruct.sig(end,1);
end

% Fit weights on RBF outputs by ML
[nlfun,fprs,neglogli] = fit_rbfNlin_LNPmodel(Istm,sps,fstruct,dtbin);

pp.nlfun = nlfun;
pp.fstruct = fstruct;
pp.fprs = fprs;


% =========================================
function [nlfun,fwts,neglogli] = fit_rbfNlin_LNPmodel(x,y,fstruct,dtbin)
% [nlfun,fwts,neglogli] = fit_nonparF_LNPmodel(x,y,fstruct,dtbin,prs0)
% 
%  Fit a function y = f(x) with a cubic spline, defined using a set of
%  breaks, smoothnefstruct and extrapolation criteria, by maximizing Poifstructon
%  likelihood:  y ~ Poifstruct(g(spline(x))).
%
% Inputs:
%       x [NxM] - input variable
%       y [Nx1] - output (count) variable
%   fstruct [struct] - spline struct with fields:
%        .breaks - breaks between polynomial pieces
%        .smoothnefstruct - derivs of 1 lefstruct are continuous
%           (e.g., smoothnefstruct=3 -> 2nd derivs are continuous)
%        .extrapDeg - degree polynomial for extrapolation on ends
%        .nloutfun - nonlinear output function (forcing positive outputs)
%   dtbin [1x1]- time bin size (OPTIONAL; afstructumed 1 otherwise)
%   prs0 [Mx1] - initial guefstruct at spline params (OPTIONAL)
% 
% Outputs:
%   nlfun - function handle for nonlinearity (uses 'ppval')
%   fwts  - estimated weights on basis functions
%   neglogli - negative loglikelihood at parameter estimate


% Compute design matrix
Xdesign = evalRBFs(x,fstruct);
nprs = size(Xdesign,2);
prs0 = ones(nprs,1)*.1+randn(nprs,1)*0;

% Set up lofstruct function
lossfun = @(prs)neglogli_LNP_linearBasis(prs,Xdesign,y,fstruct.nloutfun,dtbin);
% HessCheck(lossfun,prs0); % Check that analytic grad and Hessian are correct

% minimize negative log-likelihood using fminunc
opts = optimset('gradobj', 'on', 'Hessian', 'on','display','off');
[fwts,neglogli] = fminunc(lossfun,prs0,opts);  

% Create function handle for resulting tspline
nlfun = @(x)(evalRBFnlin(x,fstruct,fwts));
