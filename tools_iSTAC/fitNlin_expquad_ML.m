function pp = fitNlin_expquad_ML(stim, sps,filts,RefreshRate)
% pp = fitNlin_expquad_ML(Stim, sps,filts,RefreshRate)
%
% Fit exponentiated-quadratic nonlinearity associated with the iSTAC
% estimator using direct maximum-likelihood estimation of the quadratic
% parameters. 
%
% INPUTS:
%          stim [TxM] - stimulus
%          sps  [Tx1] - column vector of spike counts in each bin
%          filt [NM x nfilts] - stimulus filters
%   RefreshRate  [1]  - stimulus frame rate (frames / sec)
%             
%  OUTPUTS:
%          pp - param struct for LNP model (for simulation or computing log-likelihood)
 

% Check inputs
if nargin < 4
    RefreshRate = 1; % Stimulus frame rate (frames / sec)
end

[slen,nkx] = size(stim); % number of time bins and spatial elements in stimulus
nfilts = size(filts,2); % number of total filters
nkt = size(filts,1)/nkx; % number of time bins in filter

% -- Compute filtered resp to signal ------------------
if nfilts==1
    % single filter
    Istm = sameconv(stim,reshape(filts,nkt,nkx)); % filter stimulus with k
else
    % multiple filters
    Istm = zeros(slen,nfilts);
    for j = 1:nfilts
        Istm(:,j) = sameconv(stim,reshape(filts(:,j),nkt,nkx));
    end
    
end

% ---- Make design matrix --------------
[i,j] = find(triu(ones(nfilts)));  % get indices for above diagonal terms
Mquad = Istm(:,i).*Istm(:,j); % quadratic portion of design matrix
nq = size(Mquad,2); % number of columns in this part
Xdesign = [Mquad Istm];

% ---- Use glmfit to fit parameters ----
betas = glmfit(Xdesign,sps,'poisson');

% -- Insert fitted params into nonlinearity parameter struct ---
fprs.const = betas(1)+log(RefreshRate);  % constant
fprs.M = full(sparse(i,j,betas(2:1+nq),nfilts,nfilts));
fprs.b = betas(nq+2:end);

% --- insert into LNP parameter structure ---------------
pp.k = reshape(filts,nkt,nkx,nfilts);
pp.dc = 0; % DC constant included with exponentiated quadratic
pp.nlfun = @(x)expquadratic(x,fprs); % nonlinearity (exists in 'nlfuns' dir)
pp.fprs = fprs;
pp.RefreshRate = RefreshRate;
pp.mask = [];
pp.model = 'LNP';
pp.ktype = 'linear';
