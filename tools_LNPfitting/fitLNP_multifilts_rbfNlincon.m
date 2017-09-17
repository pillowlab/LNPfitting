function [pp,neglogli] = fitLNP_multifilts_rbfNlincon(pp,Stim,sps,optimArgs)
% [pp,neglogli] = fitLNP_multifilts_rbfNlincon(pp,Stim,sps,optimArgs)
%
% ML fitting of filters and RBF nonlinearity, with filters constrained to
% be unit vectors
%
%  INPUTS:
%             pp [1x1] - param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  OUTPUS:
%          ppnew [1x1] - new param struct (with estimated params)
%       neglogli [1x1] - negative log-likelihood at ML estimate
%
% updated: Feb 13, 2014 (JW Pillow)
 

% ===================================================
% Set optimization parameters 
defaultprs = {'Algorithm','sqp','GradObj','on','GradConstr','on','maxiter',4000,'maxfunevals',1e6};
%defaultprs = {'Largescale','off','maxiter',1,'maxfunevals',1e8};
if nargin > 3
    opts = optimset(defaultprs{:}, optimArgs{:});
else
    opts = optimset(defaultprs{:});
end
% ===================================================

% Set initial params 
[filtprs0,optPrs] = setupfitting_LNP(pp,Stim,sps);
optPrs.fstruct = pp.fstruct;
fprs0 = pp.fprs;
Loss = @(prs)(neglogli_LNP_multifilts_rbfNlin(prs,optPrs));  % loss function

% Remove DC component (last filter coeff)
filtprs0(end) = [];
nfiltprs = length(filtprs0);
nfilts = optPrs.nfilts;
prs0 = [filtprs0;fprs0];

% minimize negative log likelihood 
filtnrms = sum(reshape(filtprs0,[],nfilts).^2)  % check that these are all 1s when we start
confun = @(x)(norm1constraintfun(x,nfiltprs,nfilts));  % constraint function
[prs,neglogli,exitflag] = fmincon(Loss,prs0,[],[],[],[],[],[],confun,opts);
if (exitflag == 0)
    fprintf('fitLNP_multifilts_rbfNlin: max # evaluations or iterations exceeded (fminunc)\n');
end

% Put returned vals back into param structure
pp = reinsertFitPrs_LNP(pp,[prs(1:nfiltprs);0],optPrs);
fprs = prs(nfiltprs+1:end);
pp.nlfun = @(x)evalRBFnlin(x,pp.fstruct,fprs);
pp.fprs = fprs;

% %----------------------------------------------------
% % ------ Check analytic gradients -------
%  DerivCheck(Loss,prs0,opts);



