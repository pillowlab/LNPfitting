function [c,ceq,dc,dceq] = norm1constraintfun(x,nfiltprs,ncols)
% [c,ceq] = norm1constraintfun(x,nfiltprs,ncols)x 
%
% nonlinear function for implementing the constraint norm(x) = 1;
% (returns x'*x - 1)

c = [];  % inequality constraint
nx = length(x); % total input length
nk = nfiltprs/ncols; % number of params per filter

if ncols==1
    ceq = norm(x(1:nfiltprs))^2 - 1;
else
    xx = reshape(x(1:nfiltprs),[],ncols);
    ceq = sum(xx.^2,1)'-ones(ncols,1);
end

if nargout>2
    dc = []; % derivative of inequality constraint

    % Compute grad of equality constraints
    nxtra = nx-nfiltprs;
    if ncols == 1
	dceq = [2*x(1:nfiltprs); zeros(nxtra,1)];
    else
	blcks = mat2cell(2*xx,nk,ones(1,ncols));
	dceq = [blkdiag(blcks{:}); zeros(nxtra,ncols)];
    end
end
