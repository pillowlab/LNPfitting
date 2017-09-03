function v = gsorth(a, B);
%  Performs Graham-Schmidt orthogonalization
%
%  V = gsorth(B),  or V = gramschm(a, B)
%    if B is a matrix, then V will be an orthonormal matrix whose 
%    column vectors have been formed by successive GS orthogonlization
%
%  V = gramschm(a, B)
%    two args:  a is a vector, B a matrix.  Returns unit vector 
%    which has been orthogonalized to B.

if nargin == 2
    % first check that B is orthogonal
    B = orth(B);
    a = a./norm(a);
    v = gs(a, B);
else
    m = size(a,2);
    v = a(:,1)./norm(a(:,1));
    for j = 2:m
        v(:,j) = gs(a(:,j), v);
    end
    
end

function vnew = gs(v, B)
% Orthogonalizes v wrt B;  assumes that B is orthogonal

v = v./norm(v);
vnew = v-B*(B'*v);
if norm(vnew) > 1e-10
    vnew = vnew./norm(vnew);
else
    fprintf(1, '\n ERROR (gsorth):  vector is linearly dependent\n');
    fprintf(1, 'Returning non-unit vector\n');
end