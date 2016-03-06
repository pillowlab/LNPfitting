function [neglogli, dL, H] = neglogli_LNP_1filt(prs,OPTprs)
% [neglogli, dL, H] = neglogli_LNP1filt(prs,OPTprs)
%
% Compute negative log-likelihood of data under an LNP model with a single
% filter and fixed nonlinearity, as a function of filter parameters
%
% INPUTS:
%        prs [Nx1] - params: [kprs - weights for stimulus kernel; 
%                               dc - mean input]
%     OPTprs [1x1] - struct for extra fitting-related stuff 
%                    (created by "makeFittingStruct_LNP.m")
% 
% OUTPUTS:
%         neglogli [1x1] - negative log likelihood
%               dL [Nx1] - gradient with respect to prs
%                H [NxN] - Hessian
%
% Requires structure of data created in "makeFittingStruct_LNP.m"
% 
% updated: Jan 16, 2014 (JW Pillow)


etol = 1e-100;  % fudge factor for avoiding log(0)

% Extract some vals from OPTprs (Opt Prs);
nprs = length(prs);   % # of params in param vector
nktot = OPTprs.nkx*OPTprs.nkt;   % total # params for k
RefreshRate = OPTprs.RefreshRate; % assumed stimulus frame rate (# samples / s)

% Unpack LNP params;
kprs = prs(1:nktot);
dc = prs(nktot+1);
if nprs~=nktot+1
    error('problem with # of parameters passed in');
end

% Initialize likelihood, gradient and Hessian -----------
Iinj = OPTprs.MSTM*kprs + dc;
[rr,drr,ddrr] = OPTprs.nlfun(Iinj);

% Check for zero-values of rr and set them to small constant
iiz = find(rr <= etol);
rr(iiz) = etol; % Set value to small
drr(iiz) = 0;   % Set derivs here to 0
ddrr(iiz) = 0; 

rr = rr/RefreshRate;
Trm1 = sum(rr);  % non-spike term
Trm2 = -OPTprs.sps'*log(rr); % spike-related term
neglogli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    drr = drr/RefreshRate; % scale by the RefreshRate
    
    % ------- Non-spiking term -----------
    dLdk0 = OPTprs.MSTM'*drr;  % filter term
    dLdb0 = sum(drr); % constant term
    
    % -------  Spiking term --------------
    frac1 = drr./rr.*OPTprs.sps;
    dLdk1 = OPTprs.MSTM'*(frac1); % filter term
    dLdb1 = sum(frac1); % constant term

    dLdk = dLdk0 - dLdk1;  % grad w.r.t. filter params
    dLdb = dLdb0 - dLdb1;  % grad w.r.t. constant (dc) term

    dL = [dLdk; dLdb]; % gradient
    
end

% ---------  Compute Hessian -----------------
if nargout > 2
    ddrr = ddrr/RefreshRate; % scale by the RefreshRate

    % ------- Non-spiking term -----------
    Hk = OPTprs.MSTM'*(bsxfun(@times,OPTprs.MSTM,ddrr)); % stim filter term
    Hb = sum(ddrr); % constant (dc) term
    Hkb = OPTprs.MSTM'*ddrr; % cross term between k and b

    % -------  Spiking term --------------
    frac2 = OPTprs.sps.*(rr.*ddrr - drr.^2)./rr.^2;
    Hk= Hk - OPTprs.MSTM'*(bsxfun(@times,OPTprs.MSTM,frac2)); % stim filter term
    Hb =  Hb-sum(frac2);  % constant (dc) term
    Hkb = Hkb - OPTprs.MSTM'*frac2; % cross term between k and b

    % Assemble the terms into matrix
    H = [[Hk Hkb]; [Hkb' Hb]];  % Hessian

end
