function [neglogli, dL, H] = neglogli_LNP_multifilts_cbfNlin(prs,OPTprs)
% [neglogli, dL, H] = neglogli_LNP_multifilts_cbfNlin(prs,OPTprs)
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
%                H [NxN] - Hessian
%
% updated: Jan 24, 2014 (JW Pillow)

etol = 1e-100;  % fudge factor for avoiding log(0)

% Extract some vals from OPTprs (Opt Prs);
nprs = length(prs);   % # of params in param vector
nktot = OPTprs.nkx*OPTprs.nkt;   % total # params for each k
nfilts = OPTprs.nfilts; % number of filters in model
nkprs = nktot*nfilts;   % number of total parameters specifying filters
nfprs = nprs-nkprs;     % number of total params for basis functions
nfperfilt = nfprs/nfilts;  % number of basis function weights per filter
RefreshRate = OPTprs.RefreshRate; % assumed stimulus frame rate (# samples / s)

% Unpack LNP params;
kprs = reshape(prs(1:nkprs),nktot,[]);
fprs = prs(nkprs+1:end);

% Initialize likelihood, gradient and Hessian -----------
xproj = OPTprs.MSTM*kprs;


if (nargout <= 1)
    % ---------  Compute neglogli only --------

    % evaluate nonlinearity
    ff = evalCBFs(xproj,OPTprs.fstruct);  % evaluate basis functions
    rr = OPTprs.fstruct.nloutfun(ff*fprs); % combine basis funcs and evaluate output function
    rr = rr./RefreshRate;
    
    % Check for near-zero values of rr
    iiz = (rr <= etol);
    rr(iiz) = etol; % Set value to small
    
    % --- Compute neglogli -----
    neglogli = sum(rr) - OPTprs.sps'*log(rr); % [non-spike term] - [spike-term]

elseif (nargout == 2)
    % ---------  Compute neglogli and Gradient -----------------

    % put fprs into more useful (blkdiag) form for computing grads
    fprsmat = mat2cell(reshape(fprs,nfperfilt,nfilts),nfperfilt,ones(1,nfilts));
    fprsmat = blkdiag(fprsmat{:});
    
    % evaluate basis functions 
    [ff,dff] = evalCBFs(xproj,OPTprs.fstruct);
    fperfilt = ff*fprsmat; 
    dfperfilt = dff*fprsmat;
    f = sum(fperfilt,2);

    % evaluate outer nonlinearity and derivatives
    [rr,drr] = OPTprs.fstruct.nloutfun(f);
    rr = rr./RefreshRate; drr = drr./RefreshRate; 
    
    % Remove near-zero values of rr
    iiz = (rr <= etol);
    rr(iiz) = etol; drr(iiz)=0; 
    
    % --- Compute neglogli -----
    neglogli = sum(rr) - OPTprs.sps'*log(rr); % [non-spike term] - [spike-term]
    
    % ---- Compute Gradient ----
    dTerm1 = OPTprs.sps./rr - 1;
    % grad wrt filter params
    drdk = bsxfun(@times,dfperfilt,drr);
    dLdk = zeros(nktot,nfilts);
    for jj = 1:nfilts
	dLdk(:,jj) = -OPTprs.MSTM'*(drdk(:,jj).*dTerm1);
    end
    % grad wrt func basis weights
    dLdf = -ff'*(drr.*dTerm1);  % deriv wrt function weights
    dL = [dLdk(:); dLdf];  % gradient
    
elseif (nargout == 3)
    % ---------  Compute neglogli, grad and Hessian --------------

    % put fprs into more useful (blkdiag) form for computing grads
    fprsmat = mat2cell(reshape(fprs,nfperfilt,nfilts),nfperfilt,ones(1,nfilts));
    fprsmat = blkdiag(fprsmat{:});

    % evaluate basis functions 
    [ff,dff,ddff] = evalCBFs(xproj,OPTprs.fstruct);
    fperfilt = ff*fprsmat; 
    dfperfilt = dff*fprsmat;
    ddfperfilt = ddff*fprsmat;
    f = sum(fperfilt,2);
    
    % evaluate outer nonlinearity and derivatives
    [rr,drr,ddrr] = OPTprs.fstruct.nloutfun(f);
    rr = rr./RefreshRate; drr = drr./RefreshRate; ddrr = ddrr./RefreshRate;
    
    % Remove near-zero values of rr
    iiz = (rr <= etol);
    rr(iiz) = etol; drr(iiz)=0; ddrr(iiz)=0;
        
    % --- Compute neglogli -----
    neglogli = sum(rr) - OPTprs.sps'*log(rr); % [non-spike term] - [spike-term]
    
    % ---- Compute Gradient ----
    dTerm1 = OPTprs.sps./rr - 1;
    % grad wrt filter params
    drdk = bsxfun(@times,dfperfilt,drr);
    dLdk = zeros(nktot,nfilts);
    for jj = 1:nfilts
	dLdk(:,jj) = -OPTprs.MSTM'*(drdk(:,jj).*dTerm1);
    end
    % grad wrt func basis weights
    dLdf = -ff'*(drr.*dTerm1);  % deriv wrt function weights
    dL = [dLdk(:); dLdf];  % gradient
    
    % ---- Compute Hessian ----
    dTerm2 = OPTprs.sps./rr.^2;
    
    % ... filter 2nd derivs ...
    %ddrddkdiag = bsxfun(@times,dfperfilt.^2,ddrr) + bsxfun(@times,ddfperfilt,drr);
    Hk = zeros(nktot,nktot,nfilts,nfilts);
    
    for jj = 1:nfilts
	for ii = 1:jj
	    if jj == ii
		ddrddk = ddrr.*dfperfilt(:,jj).^2 + drr.*ddfperfilt(:,jj);
		dkk = -ddrddk.*dTerm1 + drdk(:,jj).^2.*dTerm2;
		Hk(:,:,jj,jj) = OPTprs.MSTM'*(bsxfun(@times,OPTprs.MSTM,dkk));
	    else
		ddrddk = ddrr.*dfperfilt(:,jj).*dfperfilt(:,ii);
		dkk = -ddrddk.*dTerm1 + drdk(:,jj).*drdk(:,ii).*dTerm2;
		Hk(:,:,jj,ii) = OPTprs.MSTM'*(bsxfun(@times,OPTprs.MSTM,dkk));
		Hk(:,:,ii,jj) = Hk(:,:,jj,ii);
	    end
	end
    end
    Hk = reshape(permute(Hk,[1 3 2 4]),nktot*nfilts,[]);
    
    % ... f weights 2nd derivs ...
    dww = -ddrr.*dTerm1 + drr.^2.*dTerm2;
    Hf = ff' * bsxfun(@times,ff,dww);
	
    % ... cross-terms d^2L / dk dw ....
    Hkf = zeros(nktot*nfilts,nfprs);
    for jj = 1:nfilts

	iiw = (jj-1)*nfperfilt+1:jj*nfperfilt; % relevant basis func weights
	iik = (jj-1)*nktot+1:jj*nktot; % relevant filter indices

	% % --- Slower, but easier to read -------
	% % (dr/dk)*(dr/dw) term
	% dkw = drr.*drdk(:,jj).*dTerm2;
	% Hk_drdk_drdw = bsxfun(@times,ff,dkw);
	%
	% % d^2l/(dkdw) term
	% Hk_d2rdkdw = bsxfun(@times,ff,ddrr.*dfperfilt(:,jj).*dTerm1);
	% Hk_d2rdkdw(:,iiw) = Hk_d2rdkdw(:,iiw) + bsxfun(@times,dff(:,iiw),drr.*dTerm1);
	% Hk_ff = (Hk_drdk_drdw-Hk_d2rdkdw);
	% % ---------------------------------------

	% --- Slightly Faster ---------------------
	Hk_ff = bsxfun(@times,ff, drr.*drdk(:,jj).*dTerm2 - ddrr.*dfperfilt(:,jj).*dTerm1);
	Hk_ff(:,iiw) = Hk_ff(:,iiw) - bsxfun(@times,dff(:,iiw),drr.*dTerm1);
	% -----------------------------------------
	
	Hkf(iik,:) = OPTprs.MSTM' * Hk_ff;
    end
    
    % Assemble full Hessian
    H = [Hk Hkf; Hkf', Hf];
    
end
