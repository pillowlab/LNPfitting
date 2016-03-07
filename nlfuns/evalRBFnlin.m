function f = evalRBFnlin(x,fstruct,fwts)
% f = evalRBFnlin(x,fstruct,fwts)
% 
% Evaluate nonlinearity parametrized with radial basis functions at values in x
%
% INPUT:
%      x [N x M] - each column represents a stimulus filter output (M <= 4)
%        fstruct - structure with fields:
%               .nfuncs [1 x 1] - number of functions for each input dimension
%               .ctrs [M x nfuncs] - each row is a list of centers 
%               .sigs [M x 1] - stdev of RBF along each input dimension
%   fwts [(nfuncs^M) x 1] - weights for combining basis functions
%
% OUTPUTS:
%      f [N x 1] - net output of nonparametric nonlinearity 
%
% Updated: 31 Jan 2014 (JW Pillow)


% Parse inputs
ndms = size(x,2); % number of filter outputs passed in
nfuncs = fstruct.nfuncs;
ctrs = fstruct.ctrs;
sigs = fstruct.sig;

% Check that the RBF dimensionality matches number of filter inputs
if size(ctrs,1) ~= ndms
    error('mismatch in # of basis functions and # centers in fstruct');    
end

switch ndms
    case 1,  % 1D input
        v = exp(-.5*(bsxfun(@minus,x,ctrs)/sigs).^2)*fwts;        

    case 2,  % 2D input
        [ii,jj] = ind2sub(nfuncs*ones(1,2),1:nfuncs^2);
        v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
        v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
        v = (v1(:,ii).*v2(:,jj))*fwts;

    case 3, % 3D input
        [ii,jj,kk] = ind2sub(nfuncs*ones(1,3),1:nfuncs^3);
        v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
        v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
        v3 = exp(-.5*(bsxfun(@minus,x(:,3),ctrs(3,:))./sigs(3)).^2);
        v = (v1(:,ii).*v2(:,jj).*v3(:,kk))*fwts;

    case 4, % 4D input
        [ii,jj,kk,ll] = ind2sub(nfuncs*ones(1,4),1:nfuncs^4);
        v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
        v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
        v3 = exp(-.5*(bsxfun(@minus,x(:,3),ctrs(3,:))./sigs(3)).^2);
        v4 = exp(-.5*(bsxfun(@minus,x(:,4),ctrs(4,:))./sigs(4)).^2);
        v = (v1(:,ii).*v2(:,jj).*v3(:,kk).*v4(:,ll))*fwts;

    otherwise
        error('RBF dimensionality greater than 4 not supported');
end

f = fstruct.nloutfun(v);


