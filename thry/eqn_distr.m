function [obj,ob_pcs]=eqn_distr(prs,ind_om,Hsc,phs_v,vit)
%[obj,ob_pcs]=eqn_avg_distr(prs,ind_om,Hsc,phs_v,vit);
%Objective function to min (zero);
%prs(1)Omga=ensembleFreq, prs(2:Nc)=oxi=lenN-1 X 1 offsets
%ind_om=lenN X 1 intrinsic freq, Hsc=lenH X 1 cell H fcn,
%phs_v=lenH X 1 cell phase var; vit=lenN X 1 vector, indicates cell type
% similar to eqn_avg_distr.m BUT doesn't use average (just own H)

Omga=prs(1);
oxi=prs(2:end);

lenN=length(ind_om);
rhs_p=zeros(lenN,1);

rhs_p(1)=ind_om(1)+interp1q(phs_v{vit(1)},Hsc{vit(1)},mod(oxi(1),1));
rhs_p(2)=ind_om(2)+interp1q(phs_v{vit(2)},Hsc{vit(2)},mod(-oxi(1),1))...
    +interp1q(phs_v{vit(2)},Hsc{vit(2)},mod(oxi(2)-oxi(1),1));
for k=3:(lenN-1)
    rhs_p(k,1)=...
        ind_om(k)+interp1q(phs_v{vit(k)},Hsc{vit(k)},mod(oxi(k)-oxi(k-1),1))...
        +interp1q(phs_v{vit(k)},Hsc{vit(k)},mod(oxi(k-2)-oxi(k-1),1));
end
rhs_p(lenN)=ind_om(end)+interp1q(phs_v{vit(lenN)},Hsc{vit(lenN)},mod(oxi(lenN-2)-oxi(end),1));

ob_pcs=Omga-rhs_p;

obj=sum(ob_pcs.^2);
