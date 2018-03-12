
%script to get I.C. in terms of phase variables from full model

load H_dataHexGrd %Nz, X_l, Y_l, ind_j

sclDist=sqrt(X_l.^2+Y_l.^2);
sclDist=sclDist./max(sclDist);

gp_lim=1.05; %the Euclidean length limit for which there is no gap junction
            %coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m
[G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j); 
G_cup(G_cup~=0)=0.5; %MAKE cells have same g_param
load disi

%******************Initial Conditions********************\\
load LimC_IC_all
for k=1:length(isi_L) %set IC to correct type (ind_j==k) and place IC_cl(xx,:)
    nm_k=length(IC_cl{k})*(1-sclDist(ind_j==k)); %[l,l-1,..1]=[center..perip]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    
    eta0_hetIO(k,1)=mean(nm_k./length(IC_cl{k}));
    
    nm_k=length(IC_cl{k})*(sclDist(ind_j==k)); %[l,l-1,..1]=[perip ... center]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    eta0_hetOI(k,1)=mean(nm_k./length(IC_cl{k}));
end


%repeat but now for homC &  homP
%for C, use same ind_j for heterog to get proper 'average' over cell types
% set to IC_cl{1} b/c all Center cells
for k=1:length(isi_L)
    nm_k=length(IC_cl{1})*(1-sclDist(ind_j==k)); %[l,l-1,..1]=[center..perip]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    
    eta0_homCIO(k,1)=mean(nm_k./length(IC_cl{1}));
    
    nm_k=length(IC_cl{1})*(sclDist(ind_j==k)); %[l,l-1,..1]=[perip ... center]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    eta0_homCOI(k,1)=mean(nm_k./length(IC_cl{1}));
end

%for P, use same ind_j for heterog to get proper 'average' over cell types
% set to IC_cl{length(isi_L)} b/c all Center cells
for k=1:length(isi_L)
    nm_k=length(IC_cl{length(isi_L)})*(1-sclDist(ind_j==k)); %[l,l-1,..1]=[center..perip]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    
    eta0_homPIO(k,1)=mean(nm_k./length(IC_cl{length(isi_L)}));
    
    nm_k=length(IC_cl{length(isi_L)})*(sclDist(ind_j==k)); %[l,l-1,..1]=[perip ... center]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    eta0_homPOI(k,1)=mean(nm_k./length(IC_cl{length(isi_L)}));
end

%remove all of NaN (cells not used)
eta0_hetIO=eta0_hetIO(~isnan(eta0_hetIO));
eta0_hetOI=eta0_hetOI(~isnan(eta0_hetOI));
eta0_homCIO=eta0_homCIO(~isnan(eta0_homCIO));
eta0_homCOI=eta0_homCOI(~isnan(eta0_homCOI));
eta0_homPIO=eta0_homPIO(~isnan(eta0_homPIO));
eta0_homPOI=eta0_homPOI(~isnan(eta0_homPOI));

save thry/Phase_IC eta0*
