% using 18 types of C<->P cell, with GRADual I.C.; HOM, all C
% using nn_coupVryStrHex.m to VARY the strength of gap (stronger in P)
%---- USING a hex grid (createHexGrid.m, creates H_dataHexGrd.mat)

%=========START============;

load H_dataHexGrd %has Nz, X_l Y_l, ind_j

gp_lim=1.05; %the Euclidean length limit for which there is no gap junction
            %coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m
[G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,0,ind_j); %MAKE cells near peripheral STRONGER gap

load disi
sclDist=sqrt(X_l.^2+Y_l.^2);
sclDist=sclDist./max(sclDist);
ind_j = ones(size(ind_j)); % all Center cells

g_parm=0.0005; %gapJ strength, in microS

G_cup=G_cup*g_parm;

%******************Initial Conditions********************\\
y0=zeros(Nz,31);
load LimC_IC_all
for k=1:length(isi_L) %set IC to correct type (ind_j==k) and place IC_cl(xx,:)
    nm_k=length(IC_cl{k})*(1-sclDist(ind_j==k)); %[l,l-1,..1]=[center..perip]
    nm_k=round(nm_k); %make integer; COULD be 0, watch in future
    nm_k(nm_k==0)=1;
    y0(ind_j==k,:)=IC_cl{k}(nm_k,:); % 
end

%------ set parameters----- 
load Parms_c2p %load saved Parms
Ci=C(ind_j);
Lcelli = Lcell(ind_j);
Rcelli = Rcell(ind_j);
gfNai = gfNa(ind_j);
gfKi = gfK(ind_j);
PCaLi = PCaL(ind_j);
PCaTi = PCaT(ind_j);
gKri = gKr(ind_j);
gKsi = gKs(ind_j);
gNai = gNa(ind_j);
INaKmaxi = INaKmax(ind_j);
KNaCai = KNaCa(ind_j);


%************Matlab's ode15s function*********************\\
tspan = [0:0.005:20];      % Integration time span, second
options = odeset('RelTol', 1e-06, 'AbsTol',1e-06, 'MaxStep', 0.5, 'Vectorized', 'On');   % Set numerical accuracy options for ODE solver  
tic
[time,V_out] = ode15s(@SA_fcn_Het,tspan,y0,options,Nz,G_cup,Ci,Lcelli,Rcelli,gfNai,gfKi,PCaLi,PCaTi,gKri,gKsi,gNai,INaKmaxi,KNaCai);
toc
V_out=V_out(:,1:Nz); %only voltages

save dhexhomC_io_g2_s1 g_parm gp_lim Nz X_l Y_l V_out time G_cup ind_j 