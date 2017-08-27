function nlin = fitNlin_hist1D(stim, sps, filt, RefreshRate, nbins)
% nlin = fitNlin_hist1D(stim, sps, filt, RefreshRate, nbins)
%
% Computes a histogram-based estimate of nonlinear function in a singl-filter LNP
% model.
%
% INPUTS:
%          stim [NxM] - stimulus
%          sps  [Nx1] - column vector of spike counts in each bin
%          filt [KxM] - stimulus filter
%   RefreshRate [1x1] - assumed stimulus refresh rate (frames per second)
%         nbins [1x1] - number of bins (optional)
%
% OUTPUT:
%          nlin [nbins x 2] - [xvals, yvals]
%
% The function forms histogram-based density estimators for P(z) and
% P(z|spike), where z is the filtered stimulus.  The nonlinearity estimate
% is then given by:   f  =  P(z|spike)*P(spike) / P(z)
% 
% updated: 11 April 2012 (JW Pillow)


% Check inputs
if (nargin < 4)
    RefreshRate = 100;
    fprintf('Assumed refresh rate of %d\n',RefreshRate);
end
if (nargin < 5) % set number of bins for nonlinearity
    nbins = min(25, max(round(length(stim)/4000), 12));  % ad hoc rule for deciding # bins
elseif isempty(nbins)
    nbins = min(25, max(round(length(stim)/4000), 12));  
end

% Filter stimulus with filter
if ~isempty(filt)
    xx = sameconv(stim, filt);
else
    xx = stim;
end

% Sort on values of filtered stimulus
[xxsort, srti] = sort(xx);
spsort = sps(srti);

% Estimate nonlinearities using even quantiles for x bins
xlen = length(xx);
ptsPerBin = xlen/nbins;
binedges = (0:ptsPerBin:xlen)';
iibins = [ceil(binedges(1:nbins)+1e-6), floor(binedges(2:end))];

% Compute nonlinearity by averaging x and spike count within each bin
nlin = zeros(nbins,2);  % nonlinearity
for j = 1:nbins
    ii = iibins(j,1):iibins(j,2);
    nlin(j,1) = mean(xxsort(ii));
    nlin(j,2) = mean(spsort(ii))*RefreshRate;
end

