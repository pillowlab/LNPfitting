function gg = fitNlin_iSTAC3(gg,Stim,sps)
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
nfilts = size(gg.k,3);

% -- Compute filtered resp to signal ------------------
if nfilts==1
    
    % single filter
    Istm = sameconv(Stim,gg.k)+gg.dc; % filter stimulus with k

else
    % multiple filters
    Istm = zeros(slen,nfilts);
    for j = 1:nfilts
        Istm(:,j) = sameconv(Stim,gg.k(:,:,j)); 
    end
    
end

[i,j] = find(triu(ones(nfilts)));  % get indices for above diagonal terms

Mquad = Istm(:,i).*Istm(:,j);
nq = size(Mquad,2);
Xdesign = [Mquad Istm];

% Now use glmfit to fit 
betas = glmfit(Xdesign,sps,'poisson');
fprs.const = betas(1)+log(RefreshRate);  % constant
fprs.M = full(sparse(i,j,betas(2:1+nq),nfilts,nfilts));
fprs.b = betas(nq+2:end);

gg.nlfun = @(x)expquadratic(x,fprs);
gg.fprs = fprs;
gg.dc = 0;