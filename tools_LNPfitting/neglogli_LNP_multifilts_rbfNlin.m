function [neglogli,dL] = neglogli_LNP_multifilts_rbfNlin(prs,OPTprs)
% [neglogli, dL] = neglogli_LNP_multifilts_rbfNlin(prs,OPTprs)
%
% Compute negative log-likelihood of data under an LNP model with
% (possibly) multiple filters. Optionally computes gradient and Hessian of
% negative log-likelihood w.r.t. filter and nonlinearity parameters.
%
% INPUTS:
%        prs [Nx1] - params: [kprs - weights for stimulus kernel; 
%                             fprs - weights on nonlinear basis funcs]
%     OPTprs [1x1] - struct for extra fitting-related stuff 
%                           (created by "setupfitting_LNP.m")
% 
% OUTPUTS:
%         neglogli [1x1] - negative log likelihood
%               dL [Nx1] - gradient with respect to prs
%
% updated: Jan 24, 2014 (JW Pillow)

etol = 1e-100;  % fudge factor for avoiding log(0)

% Extract some vals from OPTprs (Opt Prs);
nprs = length(prs);   % # of params in param vector
nktot = OPTprs.nkx*OPTprs.nkt;   % total # params for each k
nfilts = OPTprs.nfilts; % number of filters in model
nkprs = nktot*nfilts;   % number of total parameters specifying filters
nfprs = nprs-nkprs;     % number of total params for basis functions
RefreshRate = OPTprs.RefreshRate; % assumed stimulus frame rate (# samples / s)

% Unpack LNP params;
kprs = reshape(prs(1:nkprs),nktot,[]);
fprs = prs(nkprs+1:end);

% Initialize likelihood, gradient and Hessian -----------
xproj = OPTprs.MSTM*kprs;


if (nargout <= 1)
    % ---------  Compute neglogli only --------

    % evaluate nonlinearity
    ff = evalRBFs(xproj,OPTprs.fstruct);  % evaluate basis functions
    rr = OPTprs.fstruct.nloutfun(ff*fprs); % combine basis funcs and evaluate output function
    rr = rr./RefreshRate;
    
    % Check for near-zero values of rr
    iiz = (rr <= etol);
    rr(iiz) = etol; % Set value to small
    
    % --- Compute neglogli -----
    neglogli = sum(rr) - OPTprs.sps'*log(rr); % [non-spike term] - [spike-term]


elseif (nargout == 2)
    % ---------  Compute neglogli and Gradient -----------------
    
    % evaluate basis functions 
    [ff,dff] = evalRBFs(xproj,OPTprs.fstruct); % basis funcs 
    f = ff*fprs;  % weighted combination of basis funcs
    
    % evaluate outer nonlinearity and derivatives
    [rr,drr] = OPTprs.fstruct.nloutfun(f);
    rr = rr./RefreshRate; 
    drr = drr./RefreshRate; 
    
    % Remove near-zero values of rr
    iiz = (rr <= etol);
    rr(iiz) = etol; drr(iiz)=0; 
    
    % --- Compute neglogli -----
    neglogli = sum(rr) - OPTprs.sps'*log(rr); % [non-spike term] - [spike-term]
    
    % ---- Compute Gradient ----
    dTerm1 = OPTprs.sps./rr - 1;

    % grad wrt filter params
    dLdk = zeros(nktot,nfilts);
    for jj = 1:nfilts
	dLdk(:,jj) = -OPTprs.MSTM'*(dTerm1.*drr.*(dff(:,:,jj)*fprs));
    end
    % grad wrt func basis weights
    dLdf = -ff'*(drr.*dTerm1);  % deriv wrt function weights
    
    % Assemble terms
    dL = [dLdk(:); dLdf];  % gradient
    
else
    error('Hessian calculation not yet implemented');
end

