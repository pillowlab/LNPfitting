function f = evalCBFnlin(x,fstruct,fwts)
% f = evalCBFnlin(x,fstruct,fwts)
% 
% Evaluate nonlinearity parametrized with cylindrical basis functions at values in x
%
% INPUT:
%      x [N x M] - each column represents a stimulus filter output
%        fstruct - structure with fields 'nfuncs', 'ctrs', 'sigs'
%                  defining the basis functions
%   fwts [M*nfuncs x 1] - weights for combining basis functions
%
% OUTPUTS:
%      f [N x 1] - net output of nonparametric nonlinearity 
%
% Updated: 31 Jan 2014 (JW Pillow)


% Parse inputs
[slen,swid] =  size(x); % size of stimulus
nfuncs = fstruct.nfuncs; % number of basis functions
ctrs = fstruct.ctrs;
sig = fstruct.sig;

% Check number of centers 
if length(ctrs)~=nfuncs
    error('mismatch in # of basis functions and # centers in fstruct');
end

% Flip x to column vector, if warranted
if (size(x,1)==1) && (length(fwts)~=swid*nfuncs)
    x = x';
    [slen,swid] =  size(x); % size of stimulus
end

% Evaluate Gaussian basis vectors
f = zeros(slen,1);
for jj = 1:swid
    f = f+exp(-.5*(bsxfun(@minus,x(:,jj),ctrs)/sig).^2)*fwts((jj-1)*nfuncs+1:jj*nfuncs);
end
f = fstruct.nloutfun(f);
