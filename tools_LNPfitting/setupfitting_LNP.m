function [prs0,OPTprs] = setupfitting_LNP(pp, Stim, sps)
%  prs0 = setupfitting_LNP(pp, Stim, sps, maxsize);
%
%  Initialize params for fitting LNP model
%
%  Inputs: pp = glm param structure
%          Stim = stimulus (time along columns, other dims along rows)
%          sps = spike vector (spike count for each stim frame)
%
%  Output: prs0 = initial parameters extracted from pp
%          OPTprs = structure with optimization params
%
% updated: Jan 16, 2014 (JW Pillow)

% ---- Set stimulus filter --------------------------------------------
if strcmp(pp.ktype, 'linear')
    OPTprs = initfit_LNP(pp,Stim,sps); % standard GLM
    ktype = 1;
else
    error('unknown filter type');
end

% ====================
% Still to be implemented: other parametrizations of stim filter
% ------------------
% elseif strcmp(pp.ktype, 'bilinear')
%     initfit_stimMatrix_GLMbi(pp,Stim); % bilinearly-parametrized stim filter
%     ktype = 2;
% elseif strcmp(pp.ktype, 'blockbilinear');
%     ktype = 3;
% end
% ====================

OPTprs.nlfun = pp.nlfun;  % set nonlinearity

% ---- Extract parameter vector -------------------------------
if (ktype == 1)  % standard GLM
    prs0 = [pp.kt(:); pp.dc];

elseif (ktype == 2) % Bilinear k 
    prs0 = [pp.kt(:); pp.kx(:); pp.dc];

elseif (ktype == 3)  % mixed rank bilinearly parametrized stim filter
    error;
    
end


 % =========================================================
 function OPTprs = initfit_LNP(pp,Stim,sps)
% initfit_LNPstimMatrix(pp,Stim)
%  
% Initialize parameters relating to stimulus design matrix 

% ---- Set up filter and stim processing params ------------------- 
nkx = size(pp.k,2);  % number stim pixels (# x params)
nkt = size(pp.ktbas,2); % # time params per stim pixel (# t params)
nfilts = size(pp.k,3);
ncols = nkx*nkt;

[slen,swid] = size(Stim);

% ---- Check size of filter and width of stimulus ----------
if (nkx ~= swid)
    error('Mismatch between stim width and kernel width');
end

% ---- Filter stimulus with spatial and temporal bases -----
OPTprs.MSTM = zeros(slen,ncols);
for i = 1:nkx
    for j = 1:nkt
        OPTprs.MSTM(:,(i-1)*nkt+j) = sameconv(Stim(:,i),pp.ktbas(:,j));
    end
end

% ----- Apply mask ----------------------------------------------
iiLi = computeMask_LNP(pp.mask,slen); % compute mask (time bins to use)
OPTprs.MSTM = OPTprs.MSTM(iiLi,:);
OPTprs.sps = sps(iiLi);
OPTprs.nlfun = pp.nlfun;

% ---- Set fields of OPTprs -------------------------------------
OPTprs.nkx = nkx;
OPTprs.nkt = nkt;
OPTprs.nfilts = nfilts;
OPTprs.ktbas = pp.ktbas;
OPTprs.slen = slen;      % Total stimulus length (course bins)
OPTprs.RefreshRate = pp.RefreshRate;
