function [pp,neglogli] = fitLNP_multifilts_cbfNlincon(pp,Stim,sps,optimArgs)
% [pp,neglogli] = fitLNP_multifilts_cbfNlin(pp,Stim,sps,optimArgs)
%
% Fits LNP model parameters (filters and nonlinearity) with filters
% constrained to be unit vectors.
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
defaultprs = {'Algorithm','sqp','GradObj','on','GradConstr','on'};
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
Loss = @(prs)(neglogli_LNP_multifilts_cbfNlin(prs,optPrs));  % loss function

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
    fprintf('fitLNP_multifilts_cbfNlin: max # evaluations or iterations exceeded (fminunc)\n');
end

% % Compute Hessian of log-likelihood to obtain posterior covariance
% if nargout > 2 
%     [neglogli,~,H] = neglogli_LNP_nfilt_nonparF(prs,optPrs);
%     [ntk,nxk,nfilts] = size(pp.k);
%     nkprs = ntk*nxk;
%     Q = cell2mat(repmat(pp.ktbas,1,nxk),ntkbas,nxkbas*ones(1,nxk));
%     B = blkdiag(Q{:});
%     B = [[B, zeros(nkprs,1)]; [zeros(1,size(B,2)) 1]];  % basis for params
%     postCov = B*(H\B');
% end


% Put returned vals back into param structure
pp = reinsertFitPrs_LNP(pp,[prs(1:nfiltprs);0],optPrs);
fprs = prs(nfiltprs+1:end);
pp.nlfun = @(x)evalCBFnlin(x,pp.fstruct,fprs);
pp.fprs = fprs;

% %----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
%  DerivCheck(Loss,prs0,opts);
%  HessCheck_Elts(@Loss_LNPfilter_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=Loss_LNPfilter_logli(prs0); toc;



