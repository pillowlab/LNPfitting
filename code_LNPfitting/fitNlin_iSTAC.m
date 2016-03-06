function gg = fitNlin_iSTAC(gg,istacFilts,DD,sprate)
% ggnew = fitNlin_iSTAC(gg,istacFilts,DD,sprate)
%
% Sets the nonlinearity for an LNP model neuron to an exponentiated
% quadratic, using iSTAC parameter fits
%
%  INPUTS: 
%             gg -  LNP model param struct (with field "kk" for filters)
%     istacFilts [N x nfilts] - each column is an iSTAC filter 
%             DD -  struct for iSTAC nonlinearity (passed back from compiSTAC.m)
%         sprate -  mean spike count per bin 
%             
%  OUTPUTS:
%          ggnew [1x1] - new param struct (with estimated params)
%
%  updated: 21 Jan, 2014 (JW Pillow)
 
% Check inputs
RefreshRate = gg.RefreshRate; % Stimulus frame rate (frames / sec)
kk = gg.k;  % filters in current model
[nkt,nkx,nfilts] = size(kk); % size and number of filters in LNP model
nkprs = nkt*nkx; % number of parameters in a single filter
kk = reshape(kk,nkprs,nfilts);
bb = istacFilts\kk; % basis in iSTAC filter space

% Compute projected Gaussian distributions
mu1 = bb'*DD.mu1; % mean, spike-triggered dist
mu0 = bb'*DD.mu0; % mean, raw dist
L1 = bb'*DD.v1*bb; % covariance, spike-triggered dist
L0 = bb'*DD.v0*bb; % covariance, raw dist
kknew = reshape(istacFilts*bb,nkt,nkx,nfilts); % filters in new basis

% Compute parameters for quadratic function
invL0 = inv(L0);
invL1 = inv(L1);
fprs.M = .5*(invL0-invL1);
fprs.b = invL1*mu1-invL0*mu0;
fprs.const = .5*(mu0'*invL0*mu0 - mu1'*invL1*mu1) + ...
    .5*(logdet(L0)-logdet(L1))+log(sprate*RefreshRate);

gg.nlfun = @(x)expquadratic(x,fprs);
gg.fprs = fprs;
gg.kt = reshape(gg.ktbas\reshape(kknew,nkt,[]),[],nkx,nfilts);
gg.k = kknew;
gg.dc = 0;