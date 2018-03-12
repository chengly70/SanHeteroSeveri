function [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,strCent,ind_j)
%[G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,strCent,ind_j);
% Function returns nearest neighbor (nn) coupling matrix
%similar to nn_coupVryStrHex.m BUT now vary strength depending on 
%cell ID (not distance)
% INPUT: spatial location (X_l,Y_l) vectors, both Nz x 1 length
%   gp_lim: the Euclidean length limit for which there is no gap junction
%  strCent determines if coupling is strongest in center (1) or peripheral
%  (0); strCent=0 is the more physiological (according to Oren & Clancy '10)
%  ind_j is an Nz x 1 vector to identify cell type from (1,2,..,18)
%coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m
% OUTPUT: returns G_cup (Nz x Nz) which are all 0's & 1's; coupling matrix

if(length(X_l)~=length(Y_l) || length(X_l)~=length(ind_j))
    disp('X_l, Y_l, ind_j need to be the same length')
    return
end
Nz=length(X_l);

cllTyp=unique(ind_j); %sorted unique cell ID (ascending order)

maxCT=max(cllTyp); %largest index
minCT=min(cllTyp); %smallest index
minG=1/15; %smallest value; 1 is the largest
slp_ln=(1-minG)/(maxCT-minCT); %positive slope, for strCent=0

crrT=0; %temp variables
cpT=0;
g_val=[];

G_cup=sparse(Nz,Nz);

for j=1:Nz
    
    crrT=ind_j(j); %current (jth) cell type
    
    Eu_dist=sqrt( (X_l-X_l(j)).^2+(Y_l-Y_l(j)).^2 );
    indic_close=Eu_dist<gp_lim;
    
    indic_close(j)=0; %remove autaptic coupling
    
    if(strCent)
        %avg between v_j & ALL v_k
        g_val=.5*(-slp_ln*(crrT-maxCT)+minG) + .5*(-slp_ln*(ind_j-maxCT)+minG);
        
        G_cup(j,:)=g_val.*indic_close;
    else
        %avg between v_j & ALL v_k
        g_val=.5*(slp_ln*(crrT-minCT)+minG) + .5*(slp_ln*(ind_j-minCT)+minG);
        
        G_cup(j,:)=g_val.*indic_close;
    end
    
end
