function gg = fitNlin_iSTAC2(gg,Stim,sps)
% ggnew = fitNlin_iSTAC(gg,istacFilts,DD,sprate)
%
% Sets the nonlinearity for an LNP model neuron to an exponentiated
% quadratic, using iSTAC parameter fits
%
%  INPUTS: 
%             gg -  LNP model param struct (with field "kk" for filters)
%     istacFilts [N x nfilts] - each column is an iSTAC filter 
%             DD -  struct for iSTAC nonlinearity (passed back from compiSTAC.m)
%         sprate -  mean spike count per bin 
%             
%  OUTPUTS:
%          ggnew [1x1] - new param struct (with estimated params)
%
%  updated: 21 Jan, 2014 (JW Pillow)
 

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

[sta,stc,rawmu,rawcov] = simpleSTC(Istm,sps,1);
sprate = sum(sps)/size(Istm,1);

% Compute projected Gaussian distributions
mu1 = sta(:);
mu0 = rawmu;
L1 = stc;
L0 = rawcov;

% Compute parameters for quadratic function
invL0 = inv(L0);
invL1 = inv(L1);
fprs.M = .5*(invL0-invL1);
fprs.b = invL1*mu1-invL0*mu0;
fprs.const = .5*(mu0'*invL0*mu0 - mu1'*invL1*mu1) + ...
    .5*(logdet(L0)-logdet(L1))+log(sprate*RefreshRate);

gg.nlfun = @(x)expquadratic(x,fprs);
gg.fprs = fprs;
gg.dc = 0;