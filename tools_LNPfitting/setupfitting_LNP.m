function [prs0,OPTprs] = setupfitting_LNP(gg, Stim, sps)
%  prs0 = setupfitting_LNP(gg, Stim, sps, maxsize);
%
%  Initialize params for fitting LNP model
%
%  Inputs: gg = glm param structure
%          Stim = stimulus (time along columns, other dims along rows)
%          sps = spike vector (spike count for each stim frame)
%
%  Output: prs0 = initial parameters extracted from gg
%          OPTprs = structure with optimization params
%
% updated: Jan 16, 2014 (JW Pillow)

% ---- Set stimulus filter --------------------------------------------
if strcmp(gg.ktype, 'linear')
    OPTprs = initfit_LNP(gg,Stim,sps); % standard GLM
    ktype = 1;
else
    error('unknown filter type');
end

% ====================
% Still to be implemented: other parametrizations of stim filter
% ------------------
% elseif strcmp(gg.ktype, 'bilinear')
%     initfit_stimMatrix_GLMbi(gg,Stim); % bilinearly-parametrized stim filter
%     ktype = 2;
% elseif strcmp(gg.ktype, 'blockbilinear');
%     ktype = 3;
% end
% ====================

OPTprs.nlfun = gg.nlfun;  % set nonlinearity

% ---- Extract parameter vector -------------------------------
if (ktype == 1)  % standard GLM
    prs0 = [gg.kt(:); gg.dc];

elseif (ktype == 2) % Bilinear k 
    prs0 = [gg.kt(:); gg.kx(:); gg.dc];

elseif (ktype == 3)  % mixed rank bilinearly parametrized stim filter
    error;
    
end


 % =========================================================
 function OPTprs = initfit_LNP(gg,Stim,sps)
% initfit_LNPstimMatrix(gg,Stim)
%  
% Initialize parameters relating to stimulus design matrix 

% ---- Set up filter and stim processing params ------------------- 
nkx = size(gg.k,2);  % number stim pixels (# x params)
nkt = size(gg.ktbas,2); % # time params per stim pixel (# t params)
nfilts = size(gg.k,3);
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
        OPTprs.MSTM(:,(i-1)*nkt+j) = sameconv(Stim(:,i),gg.ktbas(:,j));
    end
end

% ----- Apply mask ----------------------------------------------
iiLi = computeMask_LNP(gg.mask,slen); % compute mask (time bins to use)
OPTprs.MSTM = OPTprs.MSTM(iiLi,:);
OPTprs.sps = sps(iiLi);
OPTprs.nlfun = gg.nlfun;

% ---- Set fields of OPTprs -------------------------------------
OPTprs.nkx = nkx;
OPTprs.nkt = nkt;
OPTprs.nfilts = nfilts;
OPTprs.ktbas = gg.ktbas;
OPTprs.slen = slen;      % Total stimulus length (course bins)
OPTprs.RefreshRate = gg.RefreshRate;
