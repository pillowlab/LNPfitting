function [pp,neglogli,postCov,H] = fitLNP_1filt_ML(pp,Stim,sps,optimArgs)
% FITLNP_1FILT_ML - fit LNP filter and dc term by ML under given nonlinearity
%
% [pp,neglogli,postCov,H] = fitLNP_1filt_ML(pp,Stim,sps,optimArgs)
% 
% Calls fminunc.  Uses gradient and Hessian to speed search. 
%
%  INPUTS: 
%             pp [1x1] -  param struct
%           Stim [NxM] - stimulus
%           sps  [Nx1] - spike count vector 
%     optimArgs [cell] - cell array of optimization params (optional),
%                        e.g., {'tolFun', '1e-12'}
%  OUTPUTS:
%          ppnew [1x1] - new param struct (with estimated params)
%       neglogli [1x1] - negative log-likelihood at ML estimate
%              H [KxK] - Hessian of log-likelihood at ML estimate (in the filter
%                        basis); first nkt elements are for stim filter) 
%
% updated: Jan 16, 2014 (JW Pillow)

% ===================================================
% Set optimization parameters 
if nargin > 3
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on');
end
% ===================================================

% Set initial params 
[prs0,optPrs] = setupfitting_LNP(pp,Stim,sps);
Loss = @(prs)(neglogli_LNP_1filt(prs,optPrs));  % loss function

% minimize negative log likelihood 
[prs,neglogli,exitflag] = fminunc(Loss,prs0,opts);
if (exitflag == 0)
    fprintf('fitLNP_1filt_ML: max # function evaluations or iterations exceeded\n');
end

% Compute Hessian of log-likelihood to obtain posterior covariance
if nargout > 2 
    [neglogli,~,H] = neglogli_LNP_1filt(prs,optPrs);
    [ntk,nxk] = size(pp.k);
    nkprs = ntk*nxk;
    % make blkdiag matrix with kt basis in each block
    [ntkbas,nxkbas] = size(pp.ktbas);
    Q = mat2cell(repmat(pp.ktbas,1,nxk),ntkbas,nxkbas*ones(1,nxk));
    B = blkdiag(Q{:}); 
    B = [[B, zeros(nkprs,1)]; [zeros(1,size(B,2)) 1]];  % basis for params
    postCov = B*(H\B');
end


% Put returned vals back into param structure
pp = reinsertFitPrs_LNP(pp,prs,optPrs);


% %----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
%  HessCheck(Loss,prs0,opts);
%  HessCheck_Elts(@Loss_LNPfilter_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=Loss_LNPfilter_logli(prs0); toc;

