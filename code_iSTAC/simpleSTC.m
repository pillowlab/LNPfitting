function [STA,STC,RawMu,RawCov] = simpleSTC(Stim, sp, nkt, CriticalSize)
%  [STA,STC,RawMu,RawCov] = simpleSTC(Stim, sp, nkt);
%
%  Computes mean and covariance of raw and spike-triggered stimuli
%
%  Input:  
%   Stim = stimulus matrix; 1st dimension = time, 2nd dimension = space
%   sp = spikes;  coloumn vector of spike count in each time bin
%                 or list of spike times
%   nkt  = number of time samples to include in the spike-triggered ensemble
%   CriticalSize - maximum number of floats to store while computing cov
%                  (smaller = slower but smaller memory requirement)
%
%  Output:
%   STA = spike-triggered average
%   STC = spike-triggered covariance (covariance around the mean);
%   RawMu = mean of raw stimulus ensemble
%   RawCov = covariance of raw ensemble
%
%  Notes:  
%   (1) doesn't compute RawMu and RawCov unless nargout > 2
%   (2) raw spike-triggered 2nd moment = STC+STA(:)*STA(:)'*nsp/(nsp-1)
%
% Copyright 2010 Pillow Lab. All rights reserved.
% $Id: simpleSTC.m 563 2011-08-04 20:42:09Z memming $

%-------- Parse inputs  ---------------------------
if nargin < 4
    CriticalSize = 1e8; % max chunk size; decrease if getting "out of memory"
end

[slen,swid] = size(Stim);  % stimulus size (time bins x spatial bins).

if (length(sp) ~= size(Stim,1))
    % Convert list of spike times to spike-count vector
    fprintf(1, 'simpleSTC: converting spike times to counts\n');
    sp = hist(sp,[1:slen])';
end
sp(1:nkt-1) = 0;  % Ignore spikes before time n


% ---------------------------------------------------
% 1. Compute only the spike-triggered STA and STC
if nargout <= 2 
    iisp = find(sp~=0);
    splen = length(iisp);
    nsp = sum(sp(iisp));
    Msz = splen*size(Stim,2)*nkt; % size of "full" spike-triggered stimulus matrix

    if Msz < CriticalSize  % Compute in one chunk if small enough
        SS = makeStimRows(Stim, nkt, iisp);
        rowlen = size(SS,2);
        STA = (sp(iisp)'*SS)'/nsp;
        if nargout > 1
            STC = SS'*(SS.*repmat(sp(iisp),1,rowlen))/(nsp-1) - STA*STA'*nsp/(nsp-1);
        end

    else % Compute in multiple chunks if too large
        nchunk = ceil(Msz/CriticalSize);
        chunksize = ceil(length(iisp)/nchunk);
        fprintf(1, 'simpleSTC: using %d chunks to compute STA/STC\n', nchunk);

        % Initialize on 1st chunk
        i0 = 1;
        imx = chunksize;
        SS = makeStimRows(Stim,nkt,iisp(i0:imx));
        rowlen = size(SS,2);
        STA = (sp(iisp(i0:imx))'*SS)';
        if nargout > 1
            STC = SS'*(SS.*repmat(sp(iisp(i0:imx)),1,rowlen));
        end
        % Compute for remaining chunks
        for j = 2:nchunk
            i0 = chunksize*(j-1)+1;
            imx = min(chunksize*j, splen);
            SS = makeStimRows(Stim,nkt,iisp(i0:imx));
            STA = STA + (sp(iisp(i0:imx))'*SS)';
            if nargout > 1
                STC = STC + SS'*(SS.*repmat(sp(iisp(i0:imx)),1,rowlen));
            end
        end
        % normalize by number of samples
        STA = STA/nsp;
        if nargout > 1
            STC = STC/(nsp-1) - STA*STA'*nsp/(nsp-1);
        end

    end
    
% -----------------------------------------------------------------------
% 2. Compute both the spike-triggered and raw stimulus means and covariances
else  

    sp = sp(nkt:end);      % Take spikes only from the nth index onward
    slen = length(sp);   % length of time indices for stim and spikes
    swid = size(Stim,2);
    nsp = sum(sp);       % number of spikes
    Msz = slen*swid*nkt;   % Size of full stimulus matrix
    rowlen = swid*nkt;     % Length of a single row of stimulus matrix

    if Msz < CriticalSize  % Check if stimulus is small enough to do in one chunk
        
        SS = makeStimRows(Stim,nkt,1);  % Convert stimulus to matrix where each row is one stim
    
        % Compute raw mean and covariance
        RawMu = mean(SS)';
        RawCov = (SS'*SS)/(slen-1)-RawMu*RawMu'*slen/(slen-1);

        % Compute spike-triggered mean and covariance
        iisp = find(sp>0);
        spvec = sp(iisp);
        STA = (spvec'*SS(iisp,:))'/nsp;
        STC = SS(iisp,:)'*(SS(iisp,:).*repmat(spvec,1,rowlen))/(nsp-1) - STA*STA'*nsp/(nsp-1);

    else  % Compute Full Stim matrix in chunks, compute mean and cov on chunks

        nchunk = ceil(Msz/CriticalSize);
        chunksize = ceil(slen/nchunk);
        fprintf(1, 'simpleSTC: using %d chunks to compute covariance\n', nchunk);

        % Compute statistics on first chunk
        SS = makeStimRows(Stim(1:chunksize+nkt-1,:),nkt,1);  % convert stimulus to "full" version
        spvec = sp(1:chunksize);
        iisp = find(spvec>0);
        RawMu = sum(SS)';
        RawCov = SS'*SS;
        STA = (spvec(iisp)'*SS(iisp,:))';
        STC = SS(iisp,:)'*(SS(iisp,:).*repmat(spvec(iisp),1,rowlen));

        % add to mean and covariance for remaining chunks
        for j = 2:nchunk;
            i0 = chunksize*(j-1)+1;  % starting index for chunk
            imax = min(slen,chunksize*j);  % ending index for chunk
            SS = makeStimRows(Stim(i0:imax+nkt-1,:),nkt,1);
            spvec = sp(i0:imax);
            iisp = find(spvec);

            RawMu = RawMu + sum(SS)';
            RawCov = RawCov +SS'*SS;
            STA = STA + (spvec(iisp)'*SS(iisp,:))';
            STC = STC + SS(iisp,:)'*(SS(iisp,:).*repmat(spvec(iisp),1,rowlen));
        end

        % divide means and covariances by number of samples
        RawMu = RawMu/slen;
        RawCov = RawCov/(slen-1) - RawMu*RawMu'*slen/(slen-1);

        STA = STA/nsp;
        STC = STC/(nsp-1) - STA*STA'*nsp/(nsp-1);
    end
end

STA = reshape(STA,[],swid);
