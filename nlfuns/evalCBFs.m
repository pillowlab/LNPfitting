function [f,df,ddf] = evalCBFs(x,fstruct)
% [f,df,ddf] = evalCBFs(x,fstruct)
% 
% Evaluate cylindrical basis functions (CBFs) with parameters specified in fstruct
%
% INPUTs: 
%    x [N x M ] - columns denote outputs of a single filter
%       fstruct - structure with basis function parameters
%                 fields: 
%                   .nfuncs - number of basis functions
%                   .ctrs - center locations
%                   .sig - stdev of CBF
%                   .includeDC - boolean for whether to include (offset) basis function
%
% OUTPUTS:
%    f [N x M*nfuncs] - each basis function, evaluated at each x
%   df [N x M*nfuncs] = 1st deriv of each basis function at each x
%  ddf [N x M*nfuncs] = 2nd deriv of each basis function at each x

% Parse inputs
[slen,swid] =  size(x);  % size of stimulus
nfuncs = fstruct.nfuncs; % number of basis functions
ncols = nfuncs*swid;     % total number of columns in f
if isfield(fstruct, 'includeDC') && (fstruct.includeDC)
    ncols = ncols+1; % increase number of columns by 1
end

ctrs = fstruct.ctrs;
sig = fstruct.sig;

if length(ctrs)~=nfuncs
    error('mismatch in # of basis functions and # centers in fstruct');
end

if nargout <=1  
    % --- evaluate Gaussian basis vectors only -----
    f = ones(slen,ncols);
    for jj = 1:swid
        f(:,(jj-1)*nfuncs+1:jj*nfuncs) = exp(-.5*(bsxfun(@minus,x(:,jj),ctrs)/sig).^2);
    end

elseif nargout == 2
    % --- evaluate Gaussian basis vectors & 1st derivs -----
    f = ones(slen,ncols);
    df = zeros(slen,ncols);
    for jj = 1:swid
        xdffs = (bsxfun(@minus,x(:,jj),ctrs));
        icols = (jj-1)*nfuncs+1:jj*nfuncs;
        f(:,icols) = exp(-.5*xdffs.^2/sig^2);
        df(:,icols) = -xdffs.*f(:,icols)/sig^2;
    end
    
elseif nargout == 3
    % --- evaluate Gaussian basis vectors & 1st and 2nd derivs -----
    f = ones(slen,ncols);
    df = zeros(slen,ncols);
    ddf = zeros(slen,ncols);
    for jj = 1:swid
        xdffs = (bsxfun(@minus,x(:,jj),ctrs));
        icols = (jj-1)*nfuncs+1:jj*nfuncs;
        f(:,icols) = exp(-.5*xdffs.^2/sig^2);
        df(:,icols) = -xdffs.*f(:,icols)/sig^2;
        ddf(:,icols) = (xdffs.^2/sig^4-1/sig^2).*f(:,icols);
    end
end