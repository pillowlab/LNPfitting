function [xgrid,ygrid,fvals,ufilts] = compNlin_2D(filts,nlfun,stim,ngridpts,fquantiles)
% [xgrid,ygrid,fvals,ufilts] = compNlin_2D(filts,nlfun,stim,ngridpts)
%
% Computes the nonlinearity in the subspace spanned by the two filters over
% the range that the stimulus occupies
%
% INPUTS:
%          filt [KxMx2] or [KM x 2] - tensor or matrix with 2 stimulus filters 
%         nlfun [func]- function handle for 2D nonlinearity
%          stim [NxM] - stimulus (optional)
%      ngridpts [1x1] - number of bins (optional: DEFAULT = 50)
%    fquantiles [1x2] - quantiles of filter responses to use for setting x
%                       and y range of plot (optional: DEFAULT = [0.025 0.975])
%
% OUTPUTS:
%      xgrid [ngrid x 1] - grid of filter 1 output values 
%      ygrid [ngrid x 1] - grid of filter 2 output values
%   vals [ngrid x ngrid] - matrix of values for nonlinearity at each grid point
%   ufilts [nfcoeff x 2] - orthonormal basis for filter subspace 


% --- Check inputs ----
if (nargin < 4) || isempty(ngridpts) % set number of bins for nonlinearity
    ngridpts = 50; % Default
end
if (nargin < 5)
    % Quantiles to use for setting range of nonlinearity to plot. 
    fquantiles = [0.025 0.975];  % Default
end

% Reshape filters and compute sizes
nxfilt = size(stim,2); % number of spatial elements in stimulus filter
if size(filts,3) == 2
   % filters are a 3-tensor
   ntfilt = size(filts,1);
   nfiltcoef = nxfilt*ntfilt; % total number of elements in each filter
else
    % filters are already an N x 2 matrix
    nfiltcoef = size(filts,1);
    ntfilt = nfiltcoef / nxfilt; 
end

% Reshape filters as vectors 
vecfilts = reshape(filts,nfiltcoef,2); %

% Convert to orthogonal unit vectors
uf1 = vecfilts(:,1)./norm(vecfilts(:,1)); % normalized filter 1
uf2 = vecfilts(:,2) - uf1*(uf1'*vecfilts(:,2)); % orthogonalize filter 2
uf2 = uf2./norm(uf2); % make unit vector of filter 2
ufilts = [uf1 uf2]; % orthonormal basis for filter subspace

% Filter stimulus with original filters to get range of normal output of filters
filtresp1 = sameconv(stim,reshape(uf1,ntfilt,nxfilt));
filtresp2 = sameconv(stim,reshape(uf2,ntfilt,nxfilt));
qf1 = quantile(filtresp1,fquantiles); % quantiles of f1 response
qf2 = quantile(filtresp2,fquantiles); % quantiles of f2 response
xgrid = linspace(qf1(1),qf1(2),ngridpts); % x grid
ygrid = linspace(qf2(1),qf2(2),ngridpts); % y grid

% Get 2D grid of points at which to evaluate nonlinearity
[xx,yy] = meshgrid(xgrid,ygrid);

% Make stimuli from basis vectors
basisStim = xx(:)*uf1' + yy(:)*uf2';

% Compute original filter responses to these basis stimuli
filtresponse = basisStim*vecfilts; 

% Evaluate nonlinearity at each point
fvals = nlfun(filtresponse); % compute nonlinearity for each point
fvals = reshape(fvals,ngridpts,ngridpts); % reshape into matrix
