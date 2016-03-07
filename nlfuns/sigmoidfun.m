function [f,df,ddf] = sigmoidfun(x,maxrate);
%  [f,df,ddf] = sigmoidfun(x,maxrate);
%
%  Implements the nonlinearity:  
%     f(x) = log(1+exp(x)).^pow;
%  Where pow = 1;
%  plus first and second derivatives

f = maxrate.*(1./(1+exp(-x)));

if nargout > 1
    df = f.*(exp(-x)./(1+exp(-x)));
end
if nargout > 2
    ddf = f.*(exp(-3*x) - exp(-x))./(1+exp(-x)).^3;
end
