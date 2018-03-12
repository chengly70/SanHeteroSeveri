

% script to plot H's from XPP; Severi et al 2012 model

g_parm=0.001; %1nS, b/c need weak coupling

load('../Parms_c2p','C') %get conductances

%---------script to plot actual H-funcs from XPP; Severi et al 2012 model
load Hfs %saved data from Hfun_[].dat

N_c=length(H);

cc=copper(N_c);


vit=[(1:2:15)';18];
figure
hold on
for j=1:length(vit)
    ind_H=vit(j);
    % must scale by C b/c of the way I calculated H in XPP (1/C on outside)
    phs_v{ind_H}=tmh{ind_H}./tmh{ind_H}(end); %so all between 0 & 1
    Hsc{ind_H}=H{ind_H}./C(ind_H)*g_parm;
    
    plot(phs_v{ind_H},Hsc{ind_H},'color',cc(ind_H,:),'LineWidth',2)
    
end
plot(tmh{1}./tmh{1}(end),zeros(size(tmh{1})),'k--')
set(gca,'FontSize',18)
set(gca,'XLim',[0 1])
xlabel('Phase')
ylabel('g*H_j')


