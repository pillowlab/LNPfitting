function [pp,lival,muLi,HLi,postCov,logEv] = fitLNP_1filt_MAP(pp,Stim,sps,CpriorInv,optimArgs)
%  [pp,lival,muLi,HLi,postCov,logEv] = fitLNP_1filt_MAP(pp,Stim,sps,CpriorInv,optimArgs)
% 
%  Computes the MAP estimate for LNP filter params, using grad and hessians.
%
%  Inputs: 
%     pp = param struct
%     Stim = stimulus
%     sps = spike vector (column vector, same height as Stim matrix)
%     CpriorInv = inverse of prior covariance matrix
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ppnew = new param struct (with ML params);
%     fval = negative log-likelihood (NOT log-posterior) at MAP estimate
%     muLi = mean of log-likelihood (based on MAP estimate)
%     HLi = inverse covariance of log-likelihood
%     postCov = posterior Covariance Matrix
%     logEv = log-evidence (at MAP value of filter)

% ===================================================
% Set optimization parameters 
if nargin > 4
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on');
end
% ===================================================

% Set initial params 
[prs0,optPrs] = setupfitting_LNP(pp,Stim,sps);
Loss = @(prs)(Loss_LNP1filt_logpost(prs,CpriorInv,optPrs));  % loss function

% minimize negative log likelihood 
[prs,~,exitflag] = fminunc(Loss,prs0,opts);
if (exitflag == 0)
    fprintf('MAPfit_LNP1D: max # function evaluations or iterations exceeded\n');
end

% Put returned vals back into param structure
pp = reinsertFitPrs_LNP(pp,prs,optPrs);

% ======= Compute (optional) additional output variables =============

% Compute Gaussian approximation to likelihood at MAP solution
if nargout > 1
    nkprs = length(prs)-1;
    [lival,~,H] = Loss_LNP1D_logli(prs,optPrs);
    HLi = H(1:nkprs,1:nkprs);
    muLi = prs(1:nkprs) + HLi\CpriorInv*prs(1:nkprs);
end

% Compute posterior covariance
if nargout > 4
    [ntk,ncols] = size(pp.k);
    nkprs = ntk*ncols;
    Q = repcell(pp.ktbas,ncols);
    B = blkdiag(Q{:});
    B = [B, zeros(nkprs,1); zeros(1,size(pp.ktbas,2)) 1];  % basis for params
    H(1:nkprs,1:nkprs) = H(1:nkprs,1:nkprs)+CpriorInv;
    postCov = B*inv(H)*B';
end

% Compute log-evidence
if nargout > 5
    kprs = prs(1:end-1);
    logEv = -lival -.5*kprs'*CpriorInv*kprs ...
        - .5*logdet(HLi*Cprior+eye(nkprs));
end


% %----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
%  HessCheck(Loss,prs0,opts);
%  HessCheck_Elts(@Loss_LNPfilter_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=Loss_LNPfilter_logli(prs0); toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  SUB-FUNCTION:  negative log-posterior (loss function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neglogpost, dL, H] = Loss_LNP1filt_logpost(prs,CpriorInv,OPTprs)
% Loss_LNP1D_logpost - negative log-posterior for LNP model filter
%
% [neglogpost, dL, H] = Loss_LNP1filt_logpost(prs,CpriorInv,OPTprs)
%
% Calls 'Loss_LNP1D_logli.m' to compute log-likelihood and adds a
% log-prior, computed from an inverse prior covariance matrix for a
% (zero-mean) Gaussian prior on the parameters.
%
% Inputs:
%          prs [Nx1] - params: [kprs - weights for stimulus kernel; 
%                               dc - mean input]
%    CpriorInv [NxN] - inverse prior covariance over params
%       OPTprs [1x1] - struct for extra fitting-related stuff 
%                     (created by setupfitting_LNP.m)
%
% Outputs:
%         neglogli [1x1] - negative log likelihood
%               dL [Nx1] - gradient with respect to prs
%                H [NxN] - Hessian
%
%   (updated: 16/03/2011 JW Pillow)


% Extract some vals from OPTprs (Opt Prs);
nktot = OPTprs.nkx*OPTprs.nkt;   % total # params for k

% Unpack LNP params;
kprs = prs(1:nktot);

switch nargout
    case 1,
        neglogli = neglogli_LNP_1filt(prs,OPTprs);
        neglogpost = neglogli + .5*kprs'*CpriorInv*kprs;
    case 2,
        [neglogli,dL] = neglogli_LNP_1filt(prs,OPTprs);
        neglogpost = neglogli + .5*kprs'*CpriorInv*kprs;
        dL = dL + CpriorInv*kprs;
    case 3
        [neglogli,dL,H] = neglogli_LNP_1filt(prs,OPTprs);
        neglogpost = neglogli + .5*kprs'*CpriorInv*kprs;
        dL(1:nktot) = dL(1:nktot) + CpriorInv*kprs;
        H(1:nktot,1:nktot) = H(1:nktot,1:nktot) + CpriorInv;
end
