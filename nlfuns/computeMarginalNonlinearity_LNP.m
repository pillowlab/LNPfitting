function [xx,rr] = computeMarginalNonlinearity_LNP(gg,filterNumber,xrnge,npts)
% computeMarginalNonlinearity_LNP - compute the effective nonlinearity
% along one filter dimension (purely for visualization purposes)
%
%  [xx,rr] = computeMarginalNonlinearity_LNP(gg,filterNumber,nstds)
%  
% Inputs:   gg [struct] - param object
%                        ------------
%                        .k - stimulus kernel
%                        .nlfun - nonlinearity
%                        ------------
%   filterNumber [1x1] - which filter to probe
%         xrange [1x1] - how far out from zero to compute nonlinearity
%
% Outputs:
%     neglogli [1x1] - negative log-likelihood of spike trains
%           rr [Nx1] - conditional intensity (in expected spikes /sec)
%      Istm [N x nK] - net linear input from stimulus (1 column per filter)
%
% updated: 21 Jan, 2014 (JW Pillow)

if nargin < 3
    xrnge = 3; % default number of standard deviations
end
if nargin < 4
    npts = 100;
end

nx = 100;  % number of points on nonlinearity
xx = linspace(-xrnge,xrnge,npts)'; % x grid for nonlinearity

nfilts = size(gg.k,3); % number of filters
kk = reshape(gg.k,[],nfilts);  % model filters
ustim = kk(:,filterNumber)./norm(kk(:,filterNumber)); % normalized stimulus in direction of filter

Stim = xx*ustim';  % stimuli spanning desired range
Istm = Stim*kk;    % filter outputs
rr = gg.nlfun(Istm);  % Conditional intensity

