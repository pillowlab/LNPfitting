function pp = reinsertFitPrs_LNP(pp,prs,OPTprs)
% pp = reinsertFitPrs_GLM(pp,prs);
%
% After fitting, reinsert params into param structure

nkt = OPTprs.nkt;  
nkx = OPTprs.nkx;   
nfilts = OPTprs.nfilts;
nktot = nkt*nkx*nfilts;

pp.kt = reshape(prs(1:nktot),nkt,nkx,nfilts);
for jj = 1:nfilts
    pp.k(:,:,jj) = pp.ktbas*pp.kt(:,:,jj);
end
pp.dc = prs(nktot+1);
