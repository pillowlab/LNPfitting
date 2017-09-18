function [pp,neglogli] = fitLNP_multifilts_rbfNlin(pp,Stim,sps,optimArgs)
% [pp,neglogli] = fitLNP_multifilts_rbfNlin(pp,Stim,sps,optimArgs)
%
%  INPUTS:
%             pp [1x1] -  param struct
%           Stim [NxM] - stimulus
%            sps [Nx1] - spike count vector 
%      optimArgs [1x1] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%
%  OUTPUS:
%          ppnew [1x1] - new param struct (with estimated params)
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
[filtprs0,optPrs] = setupfitting_LNP(pp,Stim,sps);
optPrs.fstruct = pp.fstruct;
fprs0 = pp.fprs;
Loss = @(prs)(neglogli_LNP_multifilts_rbfNlin(prs,optPrs));  % loss function

% Remove DC component (last filter coeff)
filtprs0(end) = [];
nfiltprs = length(filtprs0);
prs0 = [filtprs0;fprs0];

if ~exist('minFunc','file')
    % minimize negative log likelihood using FMINUNC
    [prs,neglogli,exitflag] = fminunc(Loss,prs0,opts);
    
    fprintf('\nMinimizing negative log-likelihood with fminunc...\n\n');
    fprintf('STRONGLY recommend installing minFunc (which is much faster)\n');
    fprintf('Get it from: https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html\n\n');
    fprintf('This code (''fitLNP_multifilts_rbfNlin.m'') will use minFunc if it exists in the path\n\n');

else
    opts.Display = 'final';
    [prs,neglogli,exitflag] = minFunc(Loss,prs0,opts);
end


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



