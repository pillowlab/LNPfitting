function pp = fitNlin_expquad_iSTACmomentbased(filts,DD,pspike,filtdims,RefreshRate)
% pp = fitNlin_expquad_iSTACmomentbased(filts,DD,pspike,filtdims,RefreshRate)
%
% Compute the exponentiated-quadratic nonlinearity for the LNP model
% nonlinearity using moments obtained directly from iSTAC estimate 
%
%  INPUTS: 
%       filts [N x nfilts] - each column is an iSTAC filter (passed back from compiSTAC.m)
%          DD - struct for iSTAC nonlinearity (passed back from compiSTAC.m)
%      pspike - mean spike count per bin 
%    filtdims - [ntdims, nxdims] specifying shape of space-time filter as matrix
% RefreshRate - stimulus frame rate (frames / sec)
%             
%  OUTPUTS:
%          pp - param struct for LNP model (for simulation or computing log-likelihood)

% Check inputs
if nargin < 4
    filtdims = [size(filts,1),1];  % assume a purely temporal filter
end
if nargin < 5
    RefreshRate = 1; % Stimulus frame rate (frames / sec)
end

% Compute parameters for quadratic function
invL0mu0 = DD.v0\DD.mu0;
invL1mu1 = DD.v1\DD.mu1;

% Compute parameters of quadratic
fprs.M = .5*(inv(DD.v0)-inv(DD.v1)); % matrix
fprs.b = invL1mu1-invL0mu0; % vector
fprs.const = .5*(DD.mu0'*invL0mu0 - DD.mu1'*invL1mu1) + ...
    .5*(logdet(DD.v0)-logdet(DD.v1))+log(pspike*RefreshRate); % constant

% insert into LNP parameter structure
nfilts = size(filts,2); % number of filters
pp.k = reshape(filts,[filtdims(:)',nfilts]);
pp.dc = 0; % DC constant included with exponentiated quadratic
pp.nlfun = @(x)expquadratic(x,fprs); % nonlinearity (exists in 'nlfuns' dir)
pp.fprs = fprs;
pp.RefreshRate = RefreshRate;
pp.mask = [];
pp.model = 'LNP';
pp.ktype = 'linear';
