function [L] = negKLsubspace(k, mu, A, bv, vA, vAv, vecs)
% [L] = negKLsubspace(k, mu, A, bv, vA,vAv); 
%  
%  loss function for computing a subapce which achieves maximal KL-divergence between a
%  Gaussian N(mu,A) and N(0,I).  
%
%  Computes KL divergence within subspace spanned by [k vecs]
%
%  inputs:  
%     k = new dimension 
%     mu = mean of Gaussian 
%     A = cov of Gaussian
%
%        Quicker to pass these in than recompute every time:
%     bV = mu' * vecs
%     vA = vecs'*A
%     vAv = vecs'*A*vecs  
%
%     vecs = basis for dimensions of subspace already "peeled off"  

if ~isempty(vecs)
    k = k - vecs*(vecs'*k); % orthogonalize k with respect to 'vecs'
end
k = k/norm(k);  %normalize k

b1 = k'*mu;
v1 = k'*A*k;
if ~isempty(bv)
    b1 = [b1; bv];
    vAb = vA*k;
    v1 = [v1 vAb'; vAb vAv];
end

L =  logdet(v1) - trace(v1) - b1'*b1;