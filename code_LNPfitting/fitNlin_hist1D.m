function [fnlin,xbinedges,xbinctrs,ybinctrs] = fitNlin_hist1D(stim, sps, filt, RefreshRate, nfbins)
% nlin = fitNlin_hist1D(stim, sps, filt, RefreshRate, nbins)
%
% Computes a histogram-based estimate of nonlinear function in a single-filter LNP model.
%
% INPUTS:
%          stim [NxM] - stimulus
%          sps  [Nx1] - column vector of spike counts in each bin
%          filt [KxM] - stimulus filter
%   RefreshRate [1x1] - assumed stimulus refresh rate (frames per second)
%        nfbins [1x1] - number of bins (optional)
%
% OUTPUT:
%          nlin - anonymous function for piecewise constant nonlinearity
%            fx - bin centers for nonlinearity
%            fy - y values for nonlinearity
%
% Compute histogram-based estimate for piece-wise constant nonlinearity by
% binning filter output and computing, for each bin, the mean spike rate

% --- Check inputs ----
if nargin < 4 || isempty(RefreshRate)
    RefreshRate = 1; 
    fprintf('fitNlin_hist1D: Assuming RefreshRate=1 (units will be spikes/bin)\n');
end
if (nargin < 5) || isempty(nfbins) % set number of bins for nonlinearity
    nfbins = min(25, max(round(length(stim)/4000), 12));  % ad hoc rule for deciding # bins
end


% Set minimum rate (so we don't produce 0 spike rate)
MINRATE = 1/length(stim);  % default choice: set to 1/number of time bins

% Filter stimulus with model filter
rawfilteroutput = sameconv(stim, filt);

% bin filter output and get bin index for each filtered stimulus
[~,xbinedges,binID] = histcounts(rawfilteroutput,nfbins); 
xbinctrs = xbinedges(1:end-1)+diff(xbinedges(1:2))/2; % use bin centers for x positions

% now compute mean spike count in each bin
ybinctrs = zeros(nfbins,1); % y values for each bin
for jj = 1:nfbins
    ybinctrs(jj) = mean(sps(binID==jj));
end
ybinctrs = max(ybinctrs,MINRATE)*RefreshRate; % divide by bin size to get units of sp/s;

% Construct function to evaluate piece-wise constant nonlinearity at any point
fnlin = @(x)(interp1(xbinctrs,ybinctrs,x,'nearest','extrap'));

