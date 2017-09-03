function gg = reinsertFitPrs_LNP(gg,prs,OPTprs)
% gg = reinsertFitPrs_GLM(gg,prs);
%
% After fitting, reinsert params into param structure

nkt = OPTprs.nkt;  
nkx = OPTprs.nkx;   
nfilts = OPTprs.nfilts;
nktot = nkt*nkx*nfilts;

gg.kt = reshape(prs(1:nktot),nkt,nkx,nfilts);
for jj = 1:nfilts
    gg.k(:,:,jj) = gg.ktbas*gg.kt(:,:,jj);
end
gg.dc = prs(nktot+1);
