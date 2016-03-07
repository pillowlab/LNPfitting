% tutorial1_MID_tfilter.m
%
% Tutorial script for illustrating MID estimator code for neuron with
% purely temporal stimulus

% 0. Set up Linear-Nonlinear-Poisson (LNP) model neuron 
nt = 32;        % number of temporal elements of filter
tvec = (-nt+1:0)'; % vector of time indices (in units of stim frames)

% Make some filters
filt1 = exp(-((tvec+nt/4)/(nt/10)).^2)-.25*exp(-((tvec+nt/2)/(nt/4)).^2); % difference of Gaussians
filt1 = filt1./norm(filt1);  %normalize

filt2 = [diff(filt1); 0];  % 2nd filter
filt2 = filt2./norm(filt2); % normalize

filt3 = [diff(filt2); 0];  % 3rd filter
filt3 = filt3./norm(filt3); % normalize

% Make plots
plot(tvec, [filt1 filt2 filt3])  
title('filters for simulation');
xlabel('time before spike'); ylabel('filter coeff');


%% 1.  Simulate data from LNP neuron

% Create stimulus 
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);
Stim = conv2(Stim,normpdf(-3:3,0,1)','same'); % smooth stimulus
RefreshRate = 100; % refresh rate

% Convolve stimulus with filters
f1 = sameconv(Stim,filt1);
f2 = sameconv(Stim,filt2);
f3 = sameconv(Stim,filt3);

% Compute output of nonlinearity
softrect = @(x)(log(1+exp(x))); % soft-rectification function
fnlin = @(x1,x2,x3)(softrect(100./(1+exp(x1-1))+10*x2.^2+4*(x3-1).^2-80));
lam = fnlin(f1,f2,f3);

%  Simulate spike train
Refreshrate = 100; % in Hz
spikes = poissrnd(lam/RefreshRate); % generate spikes


%%
[sta,stc,rawmu,rawcov] = simpleSTC(Stim,spikes,nt);
[u,s,v] = svd(stc);

ndims = 10;  % (Only need 2, but compute 10 for demonstration purposes)
eigvalthresh = 0.05; % eigenvalue cutoff threshold (for pruning dims from raw stimulus)
[vecs, vals, DD] = compiSTAC(sta, stc, rawmu, rawcov, ndims,eigvalthresh);
KLcontributed = [vals(1); diff(vals)];
ndims = length(vals);

subplot(221);  plot(1:ndims, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(221);
plot(tvec, filt1, 'k--', tvec, u(:,1:2)*u(:,1:2)'*filt1, ...
    tvec, vecs(:,1:2)*vecs(:,1:2)'*filt1, 'r');
title('Reconstruction of 1st filter'); ylabel('filter coeff');
legend('true k', 'STC', 'iSTAC', 'location', 'northwest');

subplot(223);
plot(tvec, filt2, 'k--', tvec, u(:,1:2)*u(:,1:2)'*filt2, ...
    tvec, vecs(:,1:2)*vecs(:,1:2)'*filt2, 'r');
title('Reconstruction of 2nd filter');
xlabel('time before spike'); ylabel('filter coeff');

subplot(222);  
plot(1:ndims, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(224);
plot(tvec, vecs(:,1:2));
title('iSTAC filters');
legend('1st', '2nd', 'location', 'northwest');

Errs = [subspace([filt1 filt2], u(:,1:2)) subspace([filt1 filt2], vecs(:,1:2))];
fprintf(1, 'Errors: STC=%.3f, iSTAC=%.3f\n', Errs(1), Errs(2));

