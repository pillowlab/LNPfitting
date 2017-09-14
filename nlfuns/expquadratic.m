function f = expquadratic(x,fprs)
% f = expquadratic(x,fprs)
% 
% Evaluate nonlinearity given by an exponentiated quadratic
%   f = exp( x'*M*x + b'*x + const)
% for each row of x
%
% INPUT:
%      x [N x M]  - each column represents a stimulus filter output
%   fprs [struct] - structure with fields 'M', 'b', 'const'
%
% OUTPUTS:
%      f [N x 1] - net output of nonparametric nonlinearity 

% Number of filters (inputs)
nfilts = length(fprs.b);

% Check size of input
if (size(x,2) ~= nfilts)
    % Take transpose if function is univariate 
    if (nfilts == 1) && (size(x,1) == 1)
        x = x';
    else
        error('Input doesn''t match size of quadratic form');
    end
end

% Compute function
f = exp(sum((x*fprs.M).*x,2) + x*fprs.b + fprs.const);
