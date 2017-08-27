function [neglogli,rr,Istm] = neglogli_LNP(gg,Stim,sps)
% NEGLOGLI_LNP - compute negative log-likelihood for LNP model
%
%   [neglogli,rr,Istm] = neglogli_LNP(gg,Stim,sps)
%  
% Inputs: gg [struct] - param object
%                       ------------
%                       .k - stimulus kernel
%                       .nlfun - nonlinearity
%                       ------------
%         Stim [NxM] - stimulus
%          sps [Nx1] - spikes
%
% Outputs:
%     neglogli [1x1] - negative log-likelihood of spike trains
%           rr [Nx1] - conditional intensity (in expected spikes /sec)
%      Istm [N x nK] - net linear input from stimulus (1 column per filter)
%
% updated: 16 Jan, 2014 (JW Pillow)

RefreshRate = gg.RefreshRate; % frame rate
slen = size(Stim,1); % number of time bins in stimulus

% -- Compute filtered resp to signal ------------------
if size(gg.k,3)==1
    
    % single filter
    Istm = sameconv(Stim,gg.k)+gg.dc; % filter stimulus with k

else
    % multiple filters
    nfilts = size(gg.k,3);
    Istm = zeros(slen,nfilts);
    for j = 1:nfilts
        Istm(:,j) = sameconv(Stim,gg.k(:,:,j)); 
    end
    
end

rr = gg.nlfun(Istm)/RefreshRate;  % Conditional intensity
iiLi = computeMask_LNP(gg.mask,size(rr,1));  % determine indices from mask
    
% ---- Compute negative log-likelihood ------
neglogli = -sps(iiLi)'*log(rr(iiLi)) + sum(rr(iiLi));

% If desired, pass out the conditional intensity scaled in sps/sec
if nargout > 1
    rr = rr*RefreshRate;
end

