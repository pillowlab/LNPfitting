function [logli,rr,Istm] = logli_LNP(pp,Stim,sps)
% LOGLI_LNP - compute log-likelihood for LNP model
%
%   [neglogli,rr,Istm] = logli_LNP(pp,Stim,sps)
%  
% Inputs: pp [struct] - param structure with fields:
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

RefreshRate = pp.RefreshRate; % frame rate
slen = size(Stim,1); % number of time bins in stimulus

% -- Compute filtered resp to signal ------------------
if size(pp.k,3)==1
    % single filter
    Istm = sameconv(Stim,pp.k)+pp.dc; % filter stimulus with k

else
    % multiple filters
    nfilts = size(pp.k,3);
    Istm = zeros(slen,nfilts);
    for j = 1:nfilts
        Istm(:,j) = sameconv(Stim,pp.k(:,:,j)); 
    end
    
end

% Compute conditional intensity 
rr = pp.nlfun(Istm)/RefreshRate;  % conditional intensity (spike rate) in sps/bin
iiLi = computeMask_LNP(pp.mask,size(rr,1));  % determine indices to keep from mask
    
% ---- Compute log-likelihood ------
logli = sps(iiLi)'*log(rr(iiLi)) - sum(rr(iiLi));

% If desired, pass out the conditional intensity in sps/sec
if nargout > 1
    rr = rr*RefreshRate;
end

