function d = compDklgaussian(mu1,C1,mu2,C2)
%  d = compDklgaussian(mu1, C1, mu2, C2)
%
%  Computes the KL divergence between two multivariate Gaussians
%
%  Inputs:
%    mu1 = mean of 1st Gaussian
%    C1 = covariance of 1st Gaussian
%    mu2 = mean of 2nd Gaussian
%    C2 = covariance of 2nd Gaussian
%
%  Notes:
%     D_KL = Integral(p1 log (p1/p2))
%  Analytically:  (where |*| = Determinant, Tr = trace, n = dimensionality
%     =  1/2 log(|C2|/|C1|) + 1/2 Tr(C1^.5*C2^(-1)*C1^.5)
%        + 1/2 (mu2-mu1)^T*C2^(-1)*(mu2-mu1) - 1/2 n

if nargin == 1;
    DD = mu1;
    mu1 = DD.mu1;
    mu2 = DD.mu2;
    C1 = DD.v1;
    C2 = DD.v2;
end

n = length(mu1);
b = mu2-mu1;

C1divC2 = (C2\C1);  % matrix we need

Term1 = trace(C1divC2);  % trace term
Term2 = b'*(C2\b);       % quadratic term
Term3 = -logdet(C1divC2);  % log-determinant term

d = .5*(Term1 + Term2 + Term3 - n);


