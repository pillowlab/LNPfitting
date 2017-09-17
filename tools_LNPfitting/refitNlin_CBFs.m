function [pp,neglogli] = refitNlin_CBFs(pp,Stim,sps)
% [pp,neglogli] = refitNlin_CBFs(pp,Stim,sps,fstruct)
%
% re-Fits nonlinearity in an LNP model with nonlinearity parameterized by 1D
% cylindrical basis functions (i.e., bumps for each filter output, combined
% linearly and then passed through an output nonlinearity).
%
% Differs from "fitNlin_CBFs" because it doesn't orthogonalize the filter basis or
% re-initialize the parameters for the basis functions using quantiles of the filter
% outputs 
%
%  Inputs: 
%             pp [1x1] -  param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  Outputs:
%     ppnew [1x1] - new param struct (with estimated params)
%  neglogli [1x1] - negative log-likelihood at ML estimate
%
% updated: Jan 16, 2014 (JW Pillow)

RefreshRate = pp.RefreshRate; % stimulus frame rate
[nkt,nkx,nfilts] = size(pp.k); % number of filters in LNP model 
fstruct = pp.fstruct; % has params for nonlinearity

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

% Fit spline by maximum likelihood (Newton's method)
dtbin = 1./RefreshRate;
[nlfun,fprs,neglogli] = fit_cbfNlin_LNPmodel(Istm,sps,fstruct,dtbin);

pp.nlfun = nlfun;
pp.fstruct = fstruct;
pp.fprs = fprs;


% =========================================
function [nlfun,fwts,neglogli] = fit_cbfNlin_LNPmodel(x,y,fstruct,dtbin)
% [nlfun,fwts,neglogli] = fit_cbfNlin_LNPmodel(x,y,fstruct,dtbin)
% 
%  Fit a function y = f(x) with a cubic spline, defined using a set of
%  breaks, smoothnefstruct and extrapolation criteria, by maximizing Poifstructon
%  likelihood:  y ~ Poifstruct(g(spline(x))).
%
% Inputs:
%       x [Nx1] - input variable
%       y [Nx1] - output (count) variable
%   fstruct [struct] - spline struct with fields:
%        .breaks - breaks between polynomial pieces
%        .smoothnefstruct - derivs of 1 lefstruct are continuous
%           (e.g., smoothnefstruct=3 -> 2nd derivs are continuous)
%        .extrapDeg - degree polynomial for extrapolation on ends
%        .nloutfun - nonlinear output function (forcing positive outputs)
%   dtbin [1x1]- time bin size (OPTIONAL; afstructumed 1 otherwise)
% 
% Outputs:
%   nlfun - function handle for nonlinearity (uses 'ppval')
%   fwts  - estimated weights on basis functions
%   neglogli - negative loglikelihood at parameter estimate

% Compute design matrix
Xdesign = evalCBFs(x,fstruct);
nprs = size(Xdesign,2);
prs0 = ones(nprs,1)*.1+randn(nprs,1)*.1;

% Set up lofstruct function
lossfun = @(prs)neglogli_LNP_linearBasis(prs,Xdesign,y,fstruct.nloutfun,dtbin);
% HessCheck(lossfun,prs0); % Check that analytic grad and Hessian are correct

% minimize negative log-likelihood using fminunc
opts = optimset('gradobj', 'on', 'Hessian', 'on','display','off');
[fwts,neglogli] = fminunc(lossfun,prs0,opts);  

% Create function handle for resulting tspline
nlfun = @(x)evalCBFnlin(x,fstruct,fwts);

