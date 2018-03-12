
%script to Calc avg G_cup values from the full network, gAvg_comp

%vectors to denote values to range over
svl_v=[1;2;3];    %1:P=15XC, 2:C=15XP, 3=same,.5*g_gap

%colums are g vals, calc for ALL even though won't use them all
%rows correspond to type of gap, svl_v
gAvg_comp=zeros(3,5); %first row is S1, second row is S2, 3rd row is S3


%for calc effective g for each of the 9 oscill
load ../H_dataHexGrd %has Nz, X_l Y_l, ind_j

gp_lim=1.05; %the Euclidean length limit for which there is no gap junction
            %coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m
[G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,0,ind_j); %MAKE cells near peripheral STRONGER gap

            
gvf=[0.0005;0.00125;0.002;0.003;0.004];
            
for sInd=1:3
    for gInd=1:5
        switch sInd
            case 1
                [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,0,ind_j); %MAKE cells near peripheral STRONGER gap
                G_cup=G_cup*gvf(gInd);
            case 2
                [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j); %MAKE cells near CENTER STRONGER gap
                G_cup=G_cup*gvf(gInd);
            case 3
                [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j);
                G_cup(G_cup~=0)=0.5; %MAKE cells have same g_param
                G_cup=G_cup*gvf(gInd);
        end
        
        g_nz=G_cup(G_cup~=0);
        gAvg_comp(sInd,gInd)=mean(g_nz)*1000; %so in units of nS
    end
end

save svGavgCompHex gAvg_comp