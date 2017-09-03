function [gg,neglogli] = fitLNP_multifilts_rbfNlin(gg,Stim,sps,optimArgs)
% [gg,neglogli] = fitLNP_multifilts_rbfNlin(gg,Stim,sps,optimArgs)
%
%  INPUTS:
%             gg [1x1] -  param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  OUTPUS:
%          ggnew [1x1] - new param struct (with estimated params)
%       neglogli [1x1] - negative log-likelihood at ML estimate
%
% updated: Jan 16, 2014 (JW Pillow)
 

% ===================================================
% Set optimization parameters 
defaultprs = {'Gradobj','on'};
%defaultprs = {'Largescale','off','maxiter',1,'maxfunevals',1e8};
if nargin > 3
    opts = optimset(defaultprs{:}, optimArgs{:});
else
    opts = optimset(defaultprs{:});
end
% ===================================================

% Set initial params 
[filtprs0,optPrs] = setupfitting_LNP(gg,Stim,sps);
optPrs.fstruct = gg.fstruct;
fprs0 = gg.fprs;
Loss = @(prs)(neglogli_LNP_multifilts_rbfNlin(prs,optPrs));  % loss function

% Remove DC component (last filter coeff)
filtprs0(end) = [];
nfiltprs = length(filtprs0);
prs0 = [filtprs0;fprs0];

% minimize negative log likelihood 
[prs,neglogli,exitflag] = fminunc(Loss,prs0,opts);
if (exitflag == 0)
    fprintf('fitLNP_multifilts_rbfNlin: max # evaluations or iterations exceeded (fminunc)\n');
end

% Put returned vals back into param structure
gg = reinsertFitPrs_LNP(gg,[prs(1:nfiltprs);0],optPrs);
fprs = prs(nfiltprs+1:end);
gg.nlfun = @(x)evalRBFnlin(x,gg.fstruct,fprs);
gg.fprs = fprs;

% %----------------------------------------------------
% % ------ Check analytic gradients -------
%  DerivCheck(Loss,prs0,opts);



