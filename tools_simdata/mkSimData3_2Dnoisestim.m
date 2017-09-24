% mkSimData3_2Dnoisestim.m.m
%
% Generates a simulated dataset with a response of a simulated LNP neuron
% to a 2D binary space-time stimulus. THe LNP neuron simulated has a set of
% five space-time filters, each of length 20 time bins by 16 space bins.
%
% Dataset created is used for demo script 3.
%
% Should be called from main code repository ("LNPfitting")


% Initialize paths (and create 'simdatadir' folder if necessary)
initpaths;

fname = 'simdatadir/simdata3.mat'; % pathname of file to create

%% 1. Make filters for LNP neuron

nkt = 12;  % number of time bins in filter
nkx = 16;  % number of space bins

dtBin = .01; % bin size for representing time (here .01s => 100 Hz stimulus frame rate)
RefreshRate = 1/dtBin; % Refresh rate
tk = (-nkt+1:0)'; % vector of time indices (in units of stim frames)

% Make some stimulus filters for simulated neuron
[xx,yy] = meshgrid(1:nkx,1:nkt);
fgabor = @(mx,my,sx,sy,ori,sf)(exp(-.5* ...
    (((xx-mx)/sx).^2+((yy-my)/sy).^2)) .* ... % Gaussian part
    cos((cos(ori)*(xx-mx)+sin(ori)*(yy-my))*sf/(2*pi)));

% Excitatory filters
sx = 2.5; sy = 3; ori = pi/4;  sf = 6; 
f1 = fgabor(nkx/2,nkt/2+2,sx,sy,ori,sf);   % filter 1 
f2 = fgabor(nkx/2-4,nkt/2,sx,sy,ori,sf); % filter 2 
f3 = fgabor(nkx/2+4,nkt/2,sx,sy,ori,sf); % filter 3

% Inhibitory filters
f4 = fgabor(nkx/2+3,nkt/2,sx,sy*1.2,-ori,sf);
f5 = fgabor(nkx/2-3,nkt/2,sx,sy*1.2,-ori,sf);

% % Visualize filters
% subplot(231); imagesc(f1); axis image;
% subplot(232); imagesc(f2); axis image;
% subplot(233); imagesc(f3); axis image;
% subplot(234); imagesc(f4); axis image;
% subplot(235); imagesc(f5); axis image;

% Normalize filters
f1 = f1./norm(f1(:));
f2 = f2./norm(f2(:));
f3 = f3./norm(f3(:));
f4 = f4./norm(f4(:));
f5 = f5./norm(f5(:));

filts_true = [f1(:), f2(:), f3(:), f4(:), f5(:)];

%% 2.  Simulate spike responses from LNP neuron

% Create stimulus 
slen = 2e5;  % Stimulus length (make longer to see better recovery of true filters)
Stim = (rand(slen,nkx)>0.5)*2-1; % 2D Gaussian binary white noise

% Convolve stimulus with filters
r1 = sameconv(Stim,f1); % filter 1 output
r2 = sameconv(Stim,f2); % filter 2 output 
r3 = sameconv(Stim,f3); % ...
r4 = sameconv(Stim,f4);
r5 = sameconv(Stim,f5);

% Make nonlinearity
f = @(x)(.25*log(1+exp(4*x))); % soft-rectification function
fsig = @(x)(125./(1+exp(-2*(x-2)))+1); % sigmoid rectifying at 150
fnlin = @(x1,x2,x3,x4,x5)(fsig(5*f(r1)+f(r2).^2+f(r3).^2-3*f(r4).^2-3*f(r5).^2));

% Compute output of nonlinearity
lam = fnlin(r1,r2,r3,r4,r5);

%  Simulate spike train
spikes = poissrnd(lam/RefreshRate); % generate spikes

% Inspect simulated data using STC
% --------------------------------
subplot(221); plot(lam(1:100)) % Conditional intensity
[sta,stc] = simpleSTC(Stim,spikes,nkt);
[u,s,v] = svd(stc);
subplot(222); plot(diag(s), 'o');
subplot(256); imagesc(sta);
subplot(257); imagesc(reshape(u(:,1),nkt,nkx));
subplot(258); imagesc(reshape(u(:,2),nkt,nkx));
subplot(259); imagesc(reshape(u(:,end-1),nkt,nkx));
subplot(2,5,10); imagesc(reshape(u(:,end),nkt,nkx));

% Report number of spikes and spike rate
nsp = sum(spikes); % number of spikes
fprintf('LNP simulation: %d time bins, %d spikes  (%.2f sp/s)\n', slen,nsp, nsp/slen*RefreshRate);

%% 3. Save the dataset

simdata.Stim = Stim;
simdata.spikes = spikes;
simdata.filts_true = filts_true;
simdata.dtBin = dtBin; 
simdata.RefreshRate = RefreshRate;
simdata.label = 'simulated dataset of 5-filter LNP neuron to spatio-temporal binary white noise';
save(fname, 'simdata');
