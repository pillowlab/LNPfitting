function [rbfvals,df] = evalRBFs(x,fstruct)
% [rbfvals,df] = evalRBFs(x,fstruct)
% 
% Evaluate radial basis functions at 1D, 2D or 3D points defined by the
% rows of x
%
% INPUT:
%      x [N x M] - each row is an input (M=1,2,3, or 4)
%        fstruct - structure with fields:
%               .nfuncs [1 x 1] - number of functions for each input dimension
%               .ctrs [M x nfuncs] - each row is a list of centers 
%               .sigs [M x 1] - stdev of RBF along each input dimension
%
% OUTPUTS:
%      rbfvals [N x (nfuncs^M)] - output of each of the nfuncs^M rbfs
%      df [N x (nfuncs^M) x M ] - deriv of each rbf w.r.t its M inputs
%
% Updated: 31 Jan 2014 (JW Pillow)


% Parse inputs
[slen,ndms] = size(x); % number of filter outputs passed in
nfuncs = fstruct.nfuncs;

ctrs = fstruct.ctrs;
sigs = fstruct.sig;

% Check that the RBF dimensionality matches number of filter inputs
if size(ctrs,1) ~= ndms
    error('mismatch in # of basis functions and # centers in fstruct');    
end


if nargout == 1
    
    switch ndms
	case 1,  % 1D input
	    rbfvals = exp(-.5*(bsxfun(@minus,x,ctrs)/sigs).^2);
	    
	case 2,  % 2D input
	    [ii,jj] = ind2sub(nfuncs*ones(1,2),1:nfuncs^2);
	    v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
	    v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
	    rbfvals = (v1(:,ii).*v2(:,jj));
	    
	case 3, % 3D input
	    [ii,jj,kk] = ind2sub(nfuncs*ones(1,3),1:nfuncs^3);
	    v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
	    v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
	    v3 = exp(-.5*(bsxfun(@minus,x(:,3),ctrs(3,:))./sigs(3)).^2);
	    rbfvals = (v1(:,ii).*v2(:,jj).*v3(:,kk));
	    
	case 4, % 4D input
	    [ii,jj,kk,ll] = ind2sub(nfuncs*ones(1,4),1:nfuncs^4);
	    v1 = exp(-.5*(bsxfun(@minus,x(:,1),ctrs(1,:))./sigs(1)).^2);
	    v2 = exp(-.5*(bsxfun(@minus,x(:,2),ctrs(2,:))./sigs(2)).^2);
	    v3 = exp(-.5*(bsxfun(@minus,x(:,3),ctrs(3,:))./sigs(3)).^2);
	    v4 = exp(-.5*(bsxfun(@minus,x(:,4),ctrs(4,:))./sigs(4)).^2);
	    rbfvals = (v1(:,ii).*v2(:,jj).*v3(:,kk).*v4(:,ll));
	    
	otherwise
	    error('RBF dimensionality greater than 4 not supported');
    end
    
else
    % Compute RBFs and Derivatives
    switch ndms
	case 1,  % 1D input
	    xdffs = bsxfun(@minus,x,ctrs)/sigs;
	    rbfvals = exp(-.5*xdffs.^2); % RBF values
	    df = (-xdffs.*rbfvals)/sigs;  % derivative of each
    	    
	case 2,  % 2D input
	    [ii,jj] = ind2sub(nfuncs*ones(1,2),1:nfuncs^2);
	    xdffs1 = bsxfun(@minus,x(:,1),ctrs(1,:))/sigs(1);
	    xdffs2 = bsxfun(@minus,x(:,2),ctrs(2,:))/sigs(2);
	    v1 = exp(-.5*(xdffs1.^2));
	    v2 = exp(-.5*(xdffs2.^2));
	    rbfvals = (v1(:,ii).*v2(:,jj));
	    %  Deriv
	    df = zeros(slen,nfuncs^2,ndms);
	    df(:,:,1) = -(xdffs1(:,ii).*rbfvals)/sigs(1);
	    df(:,:,2) = -(xdffs2(:,jj).*rbfvals)/sigs(2);
	    
	case 3, % 3D input
	    [ii,jj,kk] = ind2sub(nfuncs*ones(1,3),1:nfuncs^3);
	    xdffs1 = bsxfun(@minus,x(:,1),ctrs(1,:))/sigs(1);
	    xdffs2 = bsxfun(@minus,x(:,2),ctrs(2,:))/sigs(2);
	    xdffs3 = bsxfun(@minus,x(:,3),ctrs(3,:))/sigs(3);	    
	    v1 = exp(-.5*(xdffs1.^2));
	    v2 = exp(-.5*(xdffs2.^2));
	    v3 = exp(-.5*(xdffs3.^2));	    
	    rbfvals = (v1(:,ii).*v2(:,jj).*v3(:,kk));
	    %  Deriv
	    df = zeros(slen,nfuncs^3,ndms);
	    df(:,:,1) = -(xdffs1(:,ii).*rbfvals)/sigs(1);
	    df(:,:,2) = -(xdffs2(:,jj).*rbfvals)/sigs(2);
	    df(:,:,3) = -(xdffs3(:,kk).*rbfvals)/sigs(3);	    

	case 4, % 4D input
	    [ii,jj,kk,ll] = ind2sub(nfuncs*ones(1,4),1:nfuncs^4);
	    xdffs1 = bsxfun(@minus,x(:,1),ctrs(1,:))/sigs(1);
	    xdffs2 = bsxfun(@minus,x(:,2),ctrs(2,:))/sigs(2);
	    xdffs3 = bsxfun(@minus,x(:,3),ctrs(3,:))/sigs(3);	    
	    xdffs4 = bsxfun(@minus,x(:,4),ctrs(4,:))/sigs(4);	    
	    v1 = exp(-.5*(xdffs1.^2));
	    v2 = exp(-.5*(xdffs2.^2));
	    v3 = exp(-.5*(xdffs3.^2));	    
	    v4 = exp(-.5*(xdffs4.^2));
	    rbfvals = (v1(:,ii).*v2(:,jj).*v3(:,kk).*v4(:,ll));
	    %  Deriv
	    df = zeros(slen,nfuncs^4,ndms);
	    df(:,:,1) = -(xdffs1(:,ii).*rbfvals)/sigs(1);
	    df(:,:,2) = -(xdffs2(:,jj).*rbfvals)/sigs(2);
	    df(:,:,3) = -(xdffs3(:,kk).*rbfvals)/sigs(3);	    
	    df(:,:,4) = -(xdffs4(:,ll).*rbfvals)/sigs(4);	    

	    
	otherwise
	    error('RBF dimensionality greater than 4 not supported');
    end
end