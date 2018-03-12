% Augmenting the model in Severi et. al 2012 paper 
% Script to plot any one of the 18 cell types, relying on saved 
% Parameter values in load ../SA_model_network/Parms_c2p.mat,LimC_IC_all.mat,disi


jth_cell=input('Enter which type (1,2,..18): ');

%******************Initial Conditions********************\\
load LimC_IC_all.mat
[m,idC]=max(IC_cl{jth_cell}(:,1)); %start at lowest/highest volt

load disi.mat
cc=copper(18); %for plots


load Parms_c2p.mat
Cm=C;
Lclm=Lcell;
Rclm=Rcell;
gfnm=gfNa;
gfkm=gfK;
pcalm=PCaL;
pcatm=PCaT;
gkrm=gKr;
gksm=gKs;
gnam=gNa;
inakm=INaKmax;
knac=KNaCa;

%--- linearly Interp I.C. ----
y0 = IC_cl{jth_cell}(idC,:)'; 

%specify parameters FIRST so can pass into SA_fun_Interm.m
% Cell Compartments \
C = Cm(jth_cell);
Lcell = Lclm(jth_cell);
Rcell = Rclm(jth_cell);
% Sarcolemmal ion current conductance \
gfNa = gfnm(jth_cell);
gfK = gfkm(jth_cell);
PCaL = pcalm(jth_cell);
PCaT = pcatm(jth_cell);
gKr = gkrm(jth_cell);
gKs = gksm(jth_cell);
gNa = gnam(jth_cell);
INaKmax = inakm(jth_cell);
KNaCa = knac(jth_cell);


%************Matlab's ode15s function*********************\\
tspan = [0:0.001:0.8];      % Integration time span, second
options = odeset('RelTol', 1e-06, 'AbsTol',1e-06, 'MaxStep', 0.5);   % Set numerical accuracy options for ODE solver  
[time,V] = ode15s(@SA_fun_Interm,tspan,y0,options,C,Lcell,Rcell,gfNa,gfK,PCaL,PCaT,gKr,gKs,gNa,INaKmax,KNaCa);

OD=V(:,1);          % membrane voltage
ODyf=V(:,2);        % funny current gating variable
ODdL=V(:,3);        % L type Ca current activation gating variable
ODfL=V(:,4);        % L type Ca current inactivation gating variable
ODfCa=V(:,5);       % L type Ca current inactivation gating variable
ODdT=V(:,6);        % T type Ca current activation gating variable
ODfT=V(:,7);        % T type Ca current inactivation gating variable
ODpaF=V(:,8);       % IKr activation gating variable
ODpaS=V(:,9);       % IKr activation gating variable
ODPi=V(:,10);       % IKr inactivation gating variable
ODn=V(:,11);        % IKs gating variable
ODa=V(:,12);        % Ach-activated K current gating variable
ODq=V(:,13);        % Ito activation gating variable
ODr=V(:,14);        % Ito inactivativation gating variable
ODm=V(:,15);        % Na current activation gating variable 
ODh=V(:,16);        % Na current inactivation gating variable
ODR = V(:,17);      % Ca release flux from SR vis RyRs
ODO = V(:,18);      % Ca release flux from SR vis RyRs
ODI = V(:,19);      % Ca release flux from SR vis RyRs
ODRI = V(:,20);     % Ca release flux from SR vis RyRs
ODfCMi = V(:,21);   % Fractional Occupancy of calmodulin by Ca in myoplasm
ODfCMs = V(:,22);   % Fractional Occupancy of calmodulin by Ca in subspace
ODfCQ = V(:,23);    % Fractional Occupancy of calsequestrin by Ca
ODfTC = V(:,24);    % Fractional Occupancy of the troponin-Ca site by Ca
ODfTMC = V(:,25);   % Fractional Occupancy of the troponin-Mg site by Ca
ODfTMM = V(:,26);   % Fractional Occupancy of the troponin-Mg site by Mg
ODCai = V(:,27);    % Intracellular Ca concentration
ODCasub = V(:,28);  % Subspace Cca concentration
ODCansr = V(:,29);  % Ca concentration in the network SR
ODCajsr = V(:,30);  % Ca concentratioin in the junctional SR
ODNai = V(:,31);    % Intracellular Na concentration

%clear V;

%***************Calculate Current*******************************\\

%***************Cell Compartments*******************\
Lsub = 0.02;         % Distance between jSR & surface membrane(submembrane space), micrometer
Vipart = 0.46;       % Part of cell volume occupied by myoplasm
VjSRpart = 0.0012;   % Part of cell volume occupied by junctional SR
VnSRpart = 0.0116;   % Part of cell volume occupied by network SR
Vcell = 1e-9*pi*Rcell^2*Lcell;    % Cell volume, mm^3
Vsub = 1e-9*2*pi*Lsub*(Rcell-Lsub/2)*Lcell;    % Submembrane space volume, mm^3
Vi = Vipart*Vcell-Vsub;  % Myoplasmic volume, mm^3
VjSR = VjSRpart*Vcell;   % Volume of junctional SR(Ca2+ release store), mm^3
VnSR = VnSRpart*Vcell;   % Volume of network SR(Ca2+ uptake store), mm^3

%***************Fixed ion concentration (mM)********\
Ki = 140;   % intracellular K concentration
Ko = 5.4;   % extracellular K concentration 
Nao = 140;  % extracellular Na concentration
Cao = 1.8;  % extracellular Ca concentration
Mgi=2.5;    % Intracellular Mg concentration
ACh = 0;    %No Ach-activated K current

%***************Ionic Values************************\
F = 96485.3415;      % Faraday Constant, C/M
RTONF = 8314.472*310/96485.3415;       % (R*T/F) factor, mV
ENa = RTONF.*log(Nao./ODNai);          % Na equilibrium potential, mV
EK = RTONF*log(Ko/Ki);                 % K equilibrium potential, mV
ECa = 0.5.*RTONF.*log(Cao./ODCasub);   % Ca equilibrium potential, mV

%**********Sarcolemmal ion current conductance******\
gKACh = 0.00864;    % Ach-activated K current, microS 
gto = 0.002;        % Transient outward K conductance, microS 

%***Modulation of sarcolemmal ion currents by ions***\
KmKp = 1.4;         % Half-maximal Ko for INaK, mM
KmNap = 14;         % Half-maximal Nai for INaK, mM

%************NaCa Exchanger Function(mM)*************\
K1ni = 395.3;       % Intracellular Na binding to first site on NaCa  
K1no = 1628;        % Extracellular Na binding to first site on NaCa 
K2ni = 2.289;       % Intracellular Na binding to second site on NaCa
K2no = 561.4;       % Extracellular Na binding to second site on NaCa
K3ni = 26.44;       % Intracellular Na binding to third site on NaCa
K3no = 4.663;       % Extracellular Na binding to third site on NaCa
Kci = 0.0207;       % Intracellular Ca binding to NaCa transporter
Kcni = 26.44;       % Intracellular Na and Ca simultaneous binding to NaCa
Kco = 3.663;        % Extracellular Ca binding to NaCa transporter
Qci = 0.1369;       % Intracellular Ca occlusion reaction of NaCa
Qco = 0;            % Extracellular Ca occlusion reaction of NaCa
Qn = 0.4315;        % Na occlusion reaction of NaCa

%*************Ca diffusion***************************\
taodifCa = 0.00004; % Time constant of Ca diffusion from the submembrane to myoplasm, s
taotr = 0.04;       % Time constant of Ca transfer from the network to junctional SR, s

%*************SR Ca ATPase Function******************\
Kup = 0.0006;       % Half-maximal Cai for Ca uptake in the network SR, mM
Pup = 12;           % Rate constant for Ca uptake by the Ca pump in the network SR, mM/s

%*************RyR function***************************\
kiCa = 500;         % 1/(mM*S)
kim = 5;            % 1/s;
koCa = 10000;       % 1/(mM^2*s)
kom = 60;           % 1/s
ks = 25e7;          % 1/s
EC50SR = 0.45;      % mM
HSR = 2.5;          % unitless
MaxSR = 15;         % unitless
MinSR = 1;          % unitless

%*************Ca and Mg buffering*******************\
CMtot = 0.045;      % Total calmodulin concentration
CQtot = 10;         % Total calsequestrin concentration
TCtot = 0.031;      % Total concentration of the troponin-Ca2+ site
TMCtot = 0.062;     % Total concentration of the troponin-Mg2+ site
kbCM = 542;         % Ca dissociation constant for calmodulin
kbCQ = 445;         % Ca dissociation constant for calsequestrin
kbTC = 446;         % Ca dissociation constant for the troponin-Ca2+ site
kbTMC = 7.51;       % Ca dissociation constant for the troponin-Mg2+ site
kbTMM = 751;        % Mg dissociation constant for the troponin-Mg2+ site
kfCM = 227700;      % Ca association constant for calmodulin
kfCQ = 534;         % Ca association constant for calsequestrin
kfTC = 88800;       % Ca association constant for troponin
kfTMC = 227700;     % Ca association constant for troponin-Mg2+ site
kfTMM = 2277;       % Mg association constant for troponin-Mg2+ site

%*************Currents*******************************\
%funny current
Kmf = 45;           % mM
IfNa =(ODyf.^2*Ko)/(Ko+Kmf)*gfNa.*(OD-ENa);
IfK = (ODyf.^2*Ko)/(Ko+Kmf)*gfK.*(OD-EK);
If = IfNa+IfK;

%L type Ca current
IsiCa = (2.*PCaL.*OD)./(RTONF.*(1-exp((-2.*OD)./RTONF))).*(ODCasub-Cao.*exp((-2.*OD)/RTONF)).*ODdL.*ODfL.*ODfCa;
IsiK = (0.000365.*PCaL*OD)./(RTONF.*(1-exp(-OD./RTONF))).*(Ki-Ko.*exp(-OD./RTONF)).*ODdL.*ODfL.*ODfCa;
IsiNa = (0.0000185.*PCaL.*OD)./(RTONF.*(1-exp(-OD./RTONF)).*(ODNai-Nao.*exp(-OD./RTONF))).*ODdL.*ODfL.*ODfCa;
ICaL = (IsiCa+IsiK+IsiNa);

%T type Ca current
ICaT = (2.*PCaT.*OD)./(RTONF.*(1-exp(-2.*OD./RTONF))).*(ODCasub-Cao.*exp(-2.*OD./RTONF)).*ODdT.*ODfT;

%Rapidly activating delayed rectifier K current
IKr = gKr.*(OD-EK).*(0.9.*ODpaF+0.1.*ODpaS).*ODPi;

%Slowly activating delayed rectifier K current
IKs = gKs.*(OD-EK).*ODn.^2;

%Ach-Activated K current
if(ACh>0)
    IKACh = gKACh.*(OD-EK).*(1+exp((OD+20)./20)).*ODa;
else 
    IKACh = 0;
end

%Transient outward K current
Ito = gto*(OD-EK).*ODq.*ODr;

%Na current
Emh = RTONF.*log((Nao +0.12*Ko)./(ODNai +0.12*Ki)); 
INa = gNa*ODm.^3.*ODh.*(OD-Emh);

%Na-K pump current
INaK = INaKmax.*(1+(KmKp./Ko).^1.2).^(-1).*(1+(KmNap./ODNai).^1.3).^(-1).*(1+(exp(-(OD-ENa+110)./20))).^(-1);

%Na Ca Exchanger Current
di = 1+(ODCasub./Kci).*(1+exp(-Qci.*OD./RTONF)+ODNai./Kcni)+(ODNai./K1ni).*(1+(ODNai./K2ni).*(1+ODNai./K3ni));
do = 1+Cao./Kco.*(1+exp(Qco.*OD./RTONF))+Nao./K1no.*(1+Nao./K2no.*(1+Nao./K3no));

k12 = (ODCasub./Kci).*exp(-Qci.*OD./RTONF)./di;
k14 = (((ODNai./K1ni.*ODNai)./K2ni).*(1+ODNai./K3ni).*exp(Qn.*OD/(2.*RTONF)))./di;
k21 = ((Cao./Kco).*exp(Qco.*OD./RTONF))./do;
k23 = ((((Nao./K1no).*Nao)./K2no).*(1+Nao./K3no).*exp(-Qn.*OD/(2.*RTONF)))./do;
k32 = exp(Qn.*OD/(2.*RTONF));
k34 = Nao./(K3no+Nao);
k41 = exp(-Qn.*OD./(2.*RTONF));
k43 = ODNai./(K3ni+ODNai);

x1 = k41.*k34.*(k23+k21)+k21.*k32.*(k43+k41);
x2 = k32.*k43.*(k14+k12)+k41.*k12.*(k34+k32);
x3 = k14.*k43.*(k23+k21)+k12.*k23.*(k43+k41);
x4 = k23.*k34.*(k14+k12)+k14.*k21.*(k34+k32);

INaCa = (KNaCa.*(x2.*k21-x1.*k12))./(x1+x2+x3+x4);

%Total Current
Itotal = Ito+INa+If+ICaL+ICaT+IKr+IKs+IKACh+INaK+INaCa;

%*******Unit Coversion and Normalization of Currents****\
conv=1000/(C*1e6);  % Convert currents from nA to pA/pF

If=If.*conv;
ICaL=ICaL.*conv;
ICaT=ICaT.*conv;
IKr=IKr.*conv;
IKs=IKs.*conv;
IKACh=IKACh.*conv;
Ito=Ito.*conv;
INa=INa.*conv;
INaK=INaK.*conv;
INaCa=INaCa.*conv;
Itotal=Itotal.*conv;

%****************Flux************************************\
%Ca release flux from SR vis RyRs
Jrel=ks.*ODO.*(ODCajsr-ODCasub);

%Intracellular Ca fluxes
%Ca diffusion flux from submembrane space to myoplasm
Jdiff=(ODCasub-ODCai)./taodifCa;
%Ca transfer flux from the network to junctional SR
Jtr=(ODCansr-ODCajsr)./taotr;
%Ca uptake by the SR
Jup=Pup./(1+Kup./ODCai);

%***************Plots*****************************************\\

figure(39)
hold on
plot(time-isi_L(jth_cell),OD,'color',cc(jth_cell,:),'LineWidth',1);
set(gca,'FontSize',18)
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Membrane Voltage');

