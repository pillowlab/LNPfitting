function [vecs, vals, GaussParams,nullBasis] = compiSTAC(mu1,A1,mu0,A0,ndims,cthr)
% [vecs, vals, GaussParams,nullBasis] = compiSTAC(mu1,A1,mu0,A0,ndims,cthr)
%
% Computes a set of iSTAC filters -- i.e., an orthogonal basis that
% captures the maximal information about spiking given the stimuli.
% This is equal to the basis that captures maximal KL divergence between
% the two Gaussians N(mu1,A1) and N(mu0,A0).
%
% Implicitly, N(mu1,A1) is a Gaussian description of the spike-triggered
% ensemble, while N(mu0,A0) is a description of the raw stim distribution
%
% Whitens with respect to N(mu0,A0), so computation is simplified.
%
% Inputs:  
%   mu1 [n x 1] = spike-triggered average    (1/nsp*sum_i y_i*x_i)
%    A1 [n x n] = spike-triggered covariance (with mean removed)
%                          (1/nsp*sum_i y_i*x_i*x_i^T - mu1*mu1^T)
%   mu0 [n x 1] = mean of stimulus ensemble (1/N sum_i x_i)
%    A0 [n x n] = cov of stim ensemble (1/N sum_i x_i x_i - mu0*mu0^T
% ndims [1 x 1] = number of filters to estimate 
% cthr  [1 x 1] = eigenvalue threshold for whitening (OPTIONAL; DEFAULT=0.01).
%                 Will project out any dimensions for which the variance of
%                 the raw stimuli is < max(eigval)*cthr. 
%
% Ouptuts: 
%   vecs [n x ndims] = matrix with columns giving an (ordered) basis for the 
%                      maximal information-preserving subspace of degree ndims
%   vals [ndims x 1] = value of the KL divergence as subspace dimensionality increases
%                      from 1 to ndims 
%   GaussParams - structure containing the means and covariances of the two
%                 Gaussians projected into the subspace of interest.
%                (Useful if we wish to use ratio-of-Gaussians to describe the
%                 nonlinearity).
%   nullBasis [n x m] = matrix with columns giving basis for undersampled
%                       subspace of raw stimuli (i.e., which was ignored) 
%   
%
% Last updated: JW Pillow 12/2013

if nargin < 6
    cthr = 1e-2;
end

% Initialize some optimization params
vecs = [];
opts = optimset('display', 'off', 'gradobj', 'off', 'largescale', 'off', ...
    'maxfunevals', 200000, 'maxiter', 50, 'algorithm', 'Active-set');

% Compute whitening matirx
[uvecs,sdiag] = svd(A0); 
sdiag = diag(sdiag);  % eigenvalues of raw stimulus covariance
if sdiag(end)/sdiag(1) > cthr; % check condition number
    % Keep full space
    Wmat = diag(1./sqrt(sdiag))*uvecs'; % whitening matrix
    nullBasis = [];
else
    % prune some dimensions
    iikeep = sdiag>sdiag(1)*cthr;
    Wmat = diag(1./sqrt(sdiag(iikeep)))*uvecs(:,iikeep)';
    nullBasis = uvecs(:,~iikeep);
    fprintf('Pruning out %d dimensions (out of %d) from raw stimuli\n',sum(~iikeep),length(sdiag));
end
mu = Wmat*(mu1-mu0);
A =  Wmat*A1*Wmat';

% Set upper and lower bounds for optimization
nd = length(mu);
UB = ones(nd,1);
LB = -ones(nd,1);

% Compute SVD of whitened covariance, for initial guesses
[u,~] = svd(A);
a = min(ndims,floor(nd/2));
k0s = [u(:,[1:a end-a+1:end]) mu./norm(mu)];

bv = [];
vA = [];
vAv = [];

j = 1;
while j <= min(ndims,nd-1)
    BackingUP = 0;
    %  Start by finding best starting point for optimization
    kstrt = orthogonalsubset(vecs, k0s);
    v0s = 0;
    for ii = 1:size(kstrt,2);
        v0s(ii) = negKLsubspace(kstrt(:,ii),mu, A, bv, vA, vAv,vecs);
    end
    imin = find(v0s == min(v0s));  imin = imin(1);
    k0 = kstrt(:,imin);
    
    % Perform optimization -- restart if optimization doesn't terminate
    Beq = zeros(j-1,1);
    [k,~,exitflag] = fmincon(@negKLsubspace, k0,[],[],vecs',Beq,LB,UB,...
        @NormEq1,opts,mu,A,bv,vA,vAv,vecs);
    if exitflag<1  % Check convergence
        %fprintf(1, 'iSTAC-- possible error: optimization not terminated; j=%d\n',j);
        % Note: Up the optimization parameter 'niter' if worried about
        % convergence
    end
    if j > 1  % normalize k with respect to previous vecs
        k = k-vecs*(vecs'*k);
        k = k./norm(k);
    end

    % Compte KL divergence along this dimension
    vecs(:,j) = k;
    
    vals(j,1) = compDklgaussian(vecs'*mu, vecs'*A*vecs, zeros(j,1), eye(j));
    valdiffs = [vals(1); diff(vals)];
    valmarginals(j,1) = compDklgaussian(k'*mu, k'*A*k, 0, 1);
    
    % Check that vals is smaller than all previous values
    if BackingUP >= 3
        BackingUP = 0;
    elseif (valdiffs(j) > min(valdiffs(1:j-1))) & (j < nd/2)
        jj = find(valdiffs(1:j-1) < valdiffs(j));
        k0s = [k k0s];
        vecs = vecs(:,1:jj(1)-1);
        vals = vals(1:jj(1)-1);
        j = jj(1);
        fprintf(1, 'Going back to iter #%d (valdiff=%.4f)\n', j,valdiffs(end));
        BackingUP = 1;
        %
    elseif j>1 
        vv = vecs(:,[1:j-2 j]);
        valtst = compDklgaussian(vv'*mu, vv'*A*vv, zeros(j-1,1), eye(j-1));
        if valtst > vals(j-1)
            fprintf(1, 'Wrong dim possibly stripped off [%.4f %.4f]; going back to prev dim\n', ...
                vals(j-1), valtst);
            k0s = [k k0s];
            vecs = vecs(:,1:j-2);
            vals = vals(1:j-2);
            j = j-1;
            BackingUP = BackingUP+1;
        end
    end
    if ~BackingUP
        fprintf(1,' Stripped off dimension %d, KL div=[%2.4f %2.4f]\n', ...
            j, valdiffs(j), valmarginals(j));
        j = j+1;
    end

    % compute projection of A and mu onto vecs
    bv = vecs'*mu;
    vA = vecs'*A;
    vAv = vecs'*A*vecs;
end

vecs = Wmat'*vecs;
vecs = gsorth(vecs);

GaussParams.mu1 = vecs'*mu1;
GaussParams.mu0 = vecs'*mu0;
GaussParams.v1 = vecs'*A1*vecs;
GaussParams.v0 = vecs'*A0*vecs;

%  -------------------------------
function vorth = orthogonalsubset(B, vecs)
%  orthogonalize set of vectors with respect to columns of B ;
%  remove those in the column space of B

etol = 1e-10;

if isempty(B)
    vorth = vecs;
    return;
end

vorth = [];
nv = 0;
for j = 1:size(vecs,2);
    k = vecs(:,j) - B*(B\vecs(:,j));
    if norm(k) > etol;
        nv = nv+1;
        vorth(:,nv) = k./norm(k);
    end
end

%  -------------------------------
function [c,ceq,dc,dceq] = NormEq1(x, varargin);
% [c,ceq] = NormEq1(x, varargin);
%
% nonlinear function for implementing the constraint norm(x) = 1;

c = [];
ceq = x'*x-1;

if nargout>2
    dc = [];
    dceq = 2*x;
end
