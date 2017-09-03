function B = normalizecols(A)
%  NORMALIZECOLS - normalizes columns of a matrix so each is a unit vector
%
%  B = normalizecols(A);

B = bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));