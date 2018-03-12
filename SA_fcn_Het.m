function dydt = SA_fcn_Het(t,y,Nz,G_cup,C,Lcell,Rcell,gfNa,gfK,PCaL,PCaT,gKr,gKs,gNa,INaKmax,KNaCa)
%dydt = SA_fcn_Het(t,y,Nz,G_cup,
%main driver (RHS) for ode15s used in SA_netMod.m, SA_net_run.m, etc
% Allows heterogeneity in 12 parameters

%**************Stimulation Current*********************\\
Isti = 0;     

nm_vars=31; %# state vars for each cell

%********Values set to equal input values**************\\
% for j=1:Nz %loop through all Nz; save 31 state vars for each cell
% V(j,1) = y(nm_vars*(j-1)+1);           % membrane voltage
% yf(j,1) = y(nm_vars*(j-1)+2);          % funny current gating variable
% dL(j,1) = y(nm_vars*(j-1)+3);          % L type Ca current activation gating variable
% fL(j,1) = y(nm_vars*(j-1)+4);          % L type Ca current inactivation gating variable
% fCa(j,1) = y(nm_vars*(j-1)+5);         % L type Ca current inactivation gating variable
% dT(j,1) = y(nm_vars*(j-1)+6);          % T type Ca current activation gating variable
% fT(j,1) = y(nm_vars*(j-1)+7);          % T type Ca current inactivation gating variable
% paF(j,1) = y(nm_vars*(j-1)+8);         % IKr activation gating variable
% paS(j,1) = y(nm_vars*(j-1)+9);         % IKr activation gating variable
% Pi(j,1) = y(nm_vars*(j-1)+10);         % IKr inactivation gating variable
% n(j,1) = y(nm_vars*(j-1)+11);          % IKs gating variable
% a(j,1) = y(nm_vars*(j-1)+12);          % Ach-activated K current gating variable 
% q(j,1) = y(nm_vars*(j-1)+13);          % Ito activation gating variable
% r(j,1) = y(nm_vars*(j-1)+14);          % Ito inactivativation gating variable
% m(j,1) = y(nm_vars*(j-1)+15);          % Na current activation gating variable
% h(j,1) = y(nm_vars*(j-1)+16);          % Na current inactivation gating variable
% R(j,1) = y(nm_vars*(j-1)+17);          % Ca release flux from SR vis RyRs
% O(j,1) = y(nm_vars*(j-1)+18);          % Ca release flux from SR vis RyRs
% I(j,1) = y(nm_vars*(j-1)+19);          % Ca release flux from SR vis RyRs
% RI(j,1) = y(nm_vars*(j-1)+20);         % Ca release flux from SR vis RyRs
% fCMi(j,1) = y(nm_vars*(j-1)+21);       % Fractional Occupancy of calmodulin by Ca in myoplasm
% fCMs(j,1) =y(nm_vars*(j-1)+22);        % Fractional Occupancy of calmodulin by Ca in subspace
% fCQ(j,1) = y(nm_vars*(j-1)+23);        % Fractional Occupancy of calsequestrin by Ca 
% fTC(j,1) = y(nm_vars*(j-1)+24);        % Fractional Occupancy of the troponin-Ca site by Ca
% fTMC(j,1) = y(nm_vars*(j-1)+25);       % Fractional Occupancy of the troponin-Mg site by Ca
% fTMM(j,1) = y(nm_vars*(j-1)+26);       % Fractional Occupancy of the troponin-Mg site by Mg
% Cai(j,1) = y(nm_vars*(j-1)+27);        % Intracellular Ca concentration
% Casub(j,1) = y(nm_vars*(j-1)+28);      % Subspace Cca concentration
% Cansr(j,1) = y(nm_vars*(j-1)+29);      % Ca concentration in the network SR
% Cajsr(j,1) = y(nm_vars*(j-1)+30);      % Ca concentratioin in the junctional SR
% Nai(j,1) = y(nm_vars*(j-1)+31);        % Intracellular Na concentration
% end
V(1:Nz,1) = y(Nz*(1-1)+1:Nz);           % membrane voltage
yf(1:Nz,1) = y(Nz*(2-1)+1:2*Nz);          % funny current gating variable
dL(1:Nz,1) = y(Nz*(3-1)+1:3*Nz);          % L type Ca current activation gating variable
fL(1:Nz,1) = y(Nz*(4-1)+1:4*Nz);          % L type Ca current inactivation gating variable
fCa(1:Nz,1) = y(Nz*(5-1)+1:5*Nz);         % L type Ca current inactivation gating variable
dT(1:Nz,1) = y(Nz*(6-1)+1:6*Nz);          % T type Ca current activation gating variable
fT(1:Nz,1) = y(Nz*(7-1)+1:7*Nz);          % T type Ca current inactivation gating variable
paF(1:Nz,1) = y(Nz*(8-1)+1:8*Nz);         % IKr activation gating variable
paS(1:Nz,1) = y(Nz*(9-1)+1:9*Nz);         % IKr activation gating variable
Pi(1:Nz,1) = y(Nz*(10-1)+1:10*Nz);         % IKr inactivation gating variable
n(1:Nz,1) = y(Nz*(11-1)+1:11*Nz);          % IKs gating variable
a(1:Nz,1) = y(Nz*(12-1)+1:12*Nz);          % Ach-activated K current gating variable 
q(1:Nz,1) = y(Nz*(13-1)+1:13*Nz);          % Ito activation gating variable
r(1:Nz,1) = y(Nz*(14-1)+1:14*Nz);          % Ito inactivativation gating variable
m(1:Nz,1) = y(Nz*(15-1)+1:15*Nz);          % Na current activation gating variable
h(1:Nz,1) = y(Nz*(16-1)+1:16*Nz);          % Na current inactivation gating variable
R(1:Nz,1) = y(Nz*(17-1)+1:17*Nz);          % Ca release flux from SR vis RyRs
O(1:Nz,1) = y(Nz*(18-1)+1:18*Nz);          % Ca release flux from SR vis RyRs
I(1:Nz,1) = y(Nz*(19-1)+1:19*Nz);          % Ca release flux from SR vis RyRs
RI(1:Nz,1) = y(Nz*(20-1)+1:20*Nz);         % Ca release flux from SR vis RyRs
fCMi(1:Nz,1) = y(Nz*(21-1)+1:21*Nz);       % Fractional Occupancy of calmodulin by Ca in myoplasm
fCMs(1:Nz,1) =y(Nz*(22-1)+1:22*Nz);        % Fractional Occupancy of calmodulin by Ca in subspace
fCQ(1:Nz,1) = y(Nz*(23-1)+1:23*Nz);        % Fractional Occupancy of calsequestrin by Ca 
fTC(1:Nz,1) = y(Nz*(24-1)+1:24*Nz);        % Fractional Occupancy of the troponin-Ca site by Ca
fTMC(1:Nz,1) = y(Nz*(25-1)+1:25*Nz);       % Fractional Occupancy of the troponin-Mg site by Ca
fTMM(1:Nz,1) = y(Nz*(26-1)+1:26*Nz);       % Fractional Occupancy of the troponin-Mg site by Mg
Cai(1:Nz,1) = y(Nz*(27-1)+1:27*Nz);        % Intracellular Ca concentration
Casub(1:Nz,1) = y(Nz*(28-1)+1:28*Nz);      % Subspace Cca concentration
Cansr(1:Nz,1) = y(Nz*(29-1)+1:29*Nz);      % Ca concentration in the network SR
Cajsr(1:Nz,1) = y(Nz*(30-1)+1:30*Nz);      % Ca concentratioin in the junctional SR
Nai(1:Nz,1) = y(Nz*(31-1)+1:31*Nz);        % Intracellular Na concentration


%***************************Constants*************************\\
% THESE constants are identical for all Nz cells; 12 other parms are NOT
% identical
%***************Cell Compartments*******************\
%C; Membrane Capacitance, microF, passed in
%Lcell; Cell length, micrometer, passed in
Lsub = 0.02;         % Distance between jSR & surface membrane(submembrane space), micrometer
%Rcell; Cell radius, micrometer, passed in
Vipart = 0.46;       % Part of cell volume occupied by myoplasm
VjSRpart = 0.0012;   % Part of cell volume occupied by junctional SR
VnSRpart = 0.0116;   % Part of cell volume occupied by network SR
Vcell = 1e-9*pi*Rcell.^2.*Lcell;    % Cell volume, mm^3 (Nz x 1 vect)
Vsub = 1e-9*2*pi*Lsub*(Rcell-Lsub/2).*Lcell;    % Submembrane space volume, mm^3 (Nz x 1 vect)
Vi = Vipart*Vcell-Vsub;  % Myoplasmic volume, mm^3 (Nz x 1 vect)
VjSR = VjSRpart*Vcell;   % Volume of junctional SR(Ca2+ release store), mm^3 (Nz x 1 vect)
VnSR = VnSRpart*Vcell;   % Volume of network SR(Ca2+ uptake store), mm^3 (Nz x 1 vect)

%***************Fixed ion concentration (mM)********\
Ki = 140;     % intracellular K concentration
Ko = 5.4;     % extracellular K concentration 
Nao = 140;    % extracellular Na concentration
Cao = 1.8;    % extracellular Ca concentration
Mgi=2.5;      % Intracellular Mg concentration
ACh= 0.00;    % acetylcholine concentration

%***************Ionic Values************************\
F=96485.3415;     % Faraday Constant, C/M
RTONF=8314.472*310/96485.3415;    % (R*T/F) factor, mV
ENa=RTONF*log(Nao./Nai);           % Na equilibrium potential, mV VECTOR (Nz)
EK=RTONF*log(Ko/Ki);              % K equilibrium potential, mV
ECa = 0.5*RTONF*log(Cao./Casub);   % Ca equilibrium potential, mV

%**********Sarcolemmal ion current conductance******\
%gfNa; funny current Na conductance, microS, passed in
%gfK;  funny current K conductance,microS, passed in 
%PCaL; L type Ca current conductance, nA/mM, passed in
%PCaT; T type Ca current conductance, nA/mM, passed in
%gKr;  Delayed rectifier K current rapid component conductance, microS, passed in
%gKs;  Delayed rectifier K current slow component conductance, microS, passed in
gKACh = 0.00864;    % Ach-activated K current, microS 
gto = 0.002;        % Transient outward K conductance, microS 
%gNa;  Na current conductance, microS, passed in
%INaKmax; Na/K pump current, nA, passed in
%KNaCa;   Na/Ca exchanger current, nA, passed in

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
taudifCa = 0.00004; % Time constant of Ca diffusion from the submembrane to myoplasm, s
tautr = 0.04;       % Time constant of Ca transfer from the network to junctional SR, s

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
IfNa =gfNa.*(yf.^2*Ko)./(Ko+Kmf).*(V-ENa); %Nz x 1 vector
IfK = gfK.*(yf.^2*Ko)./(Ko+Kmf).*(V-EK); %Nz x 1 vector
If = IfNa+IfK; %Nz x 1 vector

%L type Ca current, all Nz x 1 vectors
IsiCa=(2*PCaL.*V)./(RTONF*(1-exp((-2*V)./RTONF))).*(Casub-Cao*exp((-2*V)./RTONF)).*dL.*fL.*fCa;
IsiK=(0.000365*PCaL.*V)./(RTONF*(1-exp(-V./RTONF))).*(Ki-Ko*exp(-V./RTONF)).*dL.*fL.*fCa;
IsiNa=(0.0000185*PCaL.*V)./(RTONF*(1-exp(-V./RTONF)).*(Nai-Nao*exp(-V./RTONF))).*dL.*fL.*fCa;
ICaL=(IsiCa+IsiK+IsiNa);

%T type Ca current, Nz x 1 vector
ICaT=(2*PCaT.*V)./(RTONF*(1-exp(-2*V/RTONF))).*(Casub-Cao*exp(-2*V./RTONF)).*dT.*fT;

%Rapidly activating delayed rectifier K current, Nz x 1 vector
IKr=gKr.*(V-EK).*(0.9*paF+0.1*paS).*Pi;

%Slowly activating delayed rectifier K current, Nz x 1 vector
IKs=gKs.*(V-EK).*n.^2;

%Ach-Activated K current, Nz x 1 vector
if(ACh>0)
    IKACh=gKACh*(V-EK).*(1+exp((V+20)./20)).*a;
else 
    IKACh=zeros(Nz,1);
end

%Transient outward K current, Nz x 1 vector
Ito=gto*(V-EK).*q.*r;

%Na current
Emh=RTONF*log((Nao +0.12*Ko)./(Nai +0.12*Ki));
INa=gNa.*m.^3.*h.*(V-Emh); %Nz x 1 vector

%Na-K pump current, Nz x 1 vector
INaK=INaKmax*(1+(KmKp/Ko)^1.2)^(-1).*(1+(KmNap./Nai).^1.3).^(-1).*(1+exp(-(V-ENa+110)./20)).^(-1);

%Na Ca Exchanger Current, Nz x 1 vector
di=1+(Casub./Kci).*(1+exp(-Qci*V./RTONF)+Nai./Kcni)+(Nai./K1ni).*(1+(Nai./K2ni).*(1+Nai./K3ni));
do=1+Cao/Kco*(1+exp(Qco*V./RTONF))+Nao/K1no*(1+Nao/K2no*(1+Nao/K3no)); 

k12=(Casub./Kci).*exp(-Qci*V./RTONF)./di;
k14=(((Nai./K1ni.*Nai)./K2ni).*(1+Nai./K3ni).*exp(Qn*V./(2*RTONF)))./di;
k21=((Cao/Kco)*exp(Qco*V./RTONF))./do;
k23=(((Nao/K1no*Nao)/K2no)*(1+Nao/K3no)*exp(-Qn*V./(2*RTONF)))./do;
k32=exp(Qn*V./(2*RTONF));
k34=Nao/(K3no+Nao); %scalar
k41=exp(-Qn*V./(2*RTONF));
k43=Nai./(K3ni+Nai);

x1=k41*k34.*(k23+k21)+k21.*k32.*(k43+k41);
x2=k32.*k43.*(k14+k12)+k41.*k12.*(k34+k32);
x3=k14.*k43.*(k23+k21)+k12.*k23.*(k43+k41);
x4=k23*k34.*(k14+k12)+k14.*k21.*(k34+k32);

INaCa=(KNaCa.*(x2.*k21-x1.*k12))./(x1+x2+x3+x4);

%****************Flux************************************\
%Ca release flux from SR vis RyRs
Jrel=ks*O.*(Cajsr-Casub);
kCaSR=MaxSR-(MaxSR-MinSR)./(1+(EC50SR./Cajsr).^HSR);
koSRCa=koCa./kCaSR;
kiSRCa=kiCa*kCaSR; %Nz x 1 vector

%Intracellular Ca fluxes
%Ca diffusion flux from submembrane space to myoplasm
Jdiff=(Casub-Cai)./taudifCa;
%Ca transfer flux from the network to junctional SR
Jtr=(Cansr-Cajsr)./tautr;
%Ca uptake by the SR
Jup=Pup./(1+Kup./Cai);

%************Ca buffering*********************************\
deltafCMi=kfCM*Cai.*(1-fCMi)-kbCM*fCMi;
deltafCMs=kfCM*Casub.*(1-fCMs)-kbCM*fCMs;
deltafCQ=kfCQ*Cajsr.*(1-fCQ)-kbCQ*fCQ;
deltafTC=kfTC*Cai.*(1-fTC)-kbTC*fTC;
deltafTMC=kfTMC*Cai.*(1-(fTMC+fTMM))-kbTMC*fTMC;
deltafTMM=kfTMM*Mgi*(1-(fTMC+fTMM))-kbTMM*fTMM;

%************SA node model Equation*********************************\\
dydt=[(1./C).*(Isti-( INa+If+Ito+ICaL+ICaT+IKr+IKs+IKACh+INaK+INaCa+sum(G_cup.*(repmat(V,1,Nz)-repmat(V',Nz,1)),2) )); % 1
        (yfss(V)-yf)./tauyf(V); 
        (dLss(V)-dL)./taudL(V);  % 3
        (fLss(V)-fL)./taufL(V); 
        (fCass(Casub)-fCa)./taufCa(Casub);  % 5                                
        (dTss(V)-dT)./taudT(V); 
        (fTss(V)-fT)./taufT(V);  % 7
        (pass(V)-paF)./taupaF(V); 
        (pass(V)-paS)./taupaS(V);  % 9
        (Piss(V)-Pi)./tauPi(V);  
        (nss(V)-n)./taun(V);  % 11
        (ass(V,ACh)-a)./taua(V,ACh);
        (qss(V)-q)./tauq(V);   % 13
        (rss(V)-r)./taur(V); 
        am(V).*(1-m)-bm(V).*m;  % 15
        ah(V).*(1-h)-bh(V).*h; 
        kim.*RI-kiSRCa.*Casub.*R-(koSRCa.*Casub.^2.*R-kom.*O);  % 17
        koSRCa.*Casub.^2.*R-kom.*O-(kiSRCa.*Casub.*O-kim.*I); 
        kiSRCa.*Casub.*O-kim.*I-(kom.*I-koSRCa.*Casub.^2.*RI);  % 19 
        kom.*I-koSRCa.*Casub.^2.*RI-(kim.*RI-kiSRCa.*Casub.*R ); 
        deltafCMi;  % 21
        deltafCMs; 
        deltafCQ;  % 23
        deltafTC; 
        deltafTMC;  % 25
        deltafTMM; 
        (Jdiff.*Vsub-Jup.*VnSR)./Vi-(CMtot.*deltafCMi+TCtot.*deltafTC+TMCtot.*deltafTMC);  % 27
        Jrel.*VjSR./Vsub-((IsiCa+ICaT-2*INaCa)./(2*F.*Vsub)+Jdiff+CMtot.*deltafCMs); 
        Jup-Jtr.*VjSR./VnSR;  % 29
        Jtr-(Jrel+CQtot.*deltafCQ); 
        -(INa+IfNa+IsiNa+3*INaK+3*INaCa)./((Vi+Vsub)*F)]; %31
% for some reason, cannot form Nz X nm_vars matrix (newline?), so use ; and reshape
% dydt=reshape(dydt,Nz,nm_vars);
% dydt=reshape(dydt',Nz*nm_vars,1);

end %end of SA_fcn_net


