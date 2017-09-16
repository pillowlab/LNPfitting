function [xgrid,fvals] = compNlin_1D(filts,nlfun,stim,ngridpts,fquantiles)
% [xgrid,fvals] = compNlin_1D(filts,nlfun,stim,ngridpts (opt),fquantiles (opt))
%
% Computes the nonlinearity in the subspace spanned by the two filters over
% the range that the stimulus occupies
%
% INPUTS:
%          filt [NxM] or [NM x 1] - tensor or matrix with 1 stimulus filter 
%         nlfun [func]- function handle for 1D nonlinearity
%          stim [NxM] - stimulus 
%      ngridpts [1x1] - number of bins (optional: DEFAULT = 50)
%    fquantiles [1x2] - quantiles of filter responses to use for setting x
%                       and y range of plot (optional: DEFAULT = [0.025 0.975])
%
% OUTPUTS:
%      xgrid [ngrid x 1] - grid of normalized filter output values 
%       vals [ngrid x 1] - matrix of values for nonlinearity at each grid point

% --- Check inputs ----
if (nargin < 4) || isempty(ngridpts) % set number of bins for nonlinearity
    ngridpts = 50; % Default
end
if (nargin < 5)
    % Quantiles to use for setting range of nonlinearity to plot. 
    fquantiles = [0.025 0.975];  % Default
end

% Reshape filters and compute sizes
nkx = size(stim,2); % number of spatial elements in stimulus filter
nkt = size(filts,1);

% Convert filter to unit vector
uf1 = filts(:)./norm(filts(:)); % normalized filter 

% Filter stimulus with original filters to get range of normal output of filters
filtresp1 = sameconv(stim,reshape(uf1,nkt,nkx));
qf1 = quantile(filtresp1,fquantiles); % quantiles of f1 response
xgrid = linspace(qf1(1),qf1(2),ngridpts); % x grid

% Make stimuli from basis vector
basisStim = xgrid(:)*uf1';

% Compute original filter responses to these basis stimuli
filtresponse = basisStim*filts(:);

% Evaluate nonlinearity at each point
fvals = nlfun(filtresponse); % compute nonlinearity for each point
