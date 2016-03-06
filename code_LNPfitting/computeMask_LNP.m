function iiLi = computeMask_LNP(mask,slen)
% iiLi = computeMask_LNP(mask,slen)
%
% Extracts spike times and bin indices for computing likelihood, given
% intervals contained in mask (an n x 2 matrix)
%
% INPUT:
%   mask [N x 2] - each row has indices [i1 i2] specifying that indices i1:i2
%                will be used for computing log-likelihood
%   slen [1 x 1] - total length of stimulus (upper bound on i2)

% Generate mask
if isempty(mask)  % No mask
    iiLi = 1:slen;
else  % Compute mask for time bins and spikes to use

    % Check that mask doesn't extend beyond stimulus window
    if max(mask(:)) > slen
        warning('GLM:mask', 'Mask is too long for data segment: truncating...');
        rowmax = max(mask,[],2);
        maskkeep = (rowmax <= slen);
        mask = mask(maskkeep,:);
    end
    
    % Construct mask
    masklens = diff(mask,1,2);
    iiLi = zeros(sum(masklens),1);
    icum = 0;
    for j = 1:size(mask,1)
        iiLi(icum+1:icum+masklens(j)) = (mask(j,1)+1:mask(j,2))';
        icum = icum+masklens(j);
    end
end