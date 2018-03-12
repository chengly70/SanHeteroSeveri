
%sim test on phase oscillators; specificy distribut of heterog

load Hfs %the uncoupled H's, must be the same size as isi_L
lenTyp=9; %using subset of all length(H)
lenH=length(H); 

cd ..; SETPATHS; cd thry/
%for calc effective g for each of the 9 oscill
load ../H_dataHexGrd %has Nz, X_l Y_l, ind_j

gp_lim=1.05; %the Euclidean length limit for which there is no gap junction
            %coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m

hetOrHom=input('Heterog (1) or Homog (2=Cent, 3=Periph)? : ');
switch hetOrHom
    case 1
        n_byt=[repmat([1;0],8,1); 0; 1]; %het type: 1,3,5,..,15,18
    case 2
        n_byt=[lenTyp;zeros(lenH-1,1)];  %all Cent
    case 3
        n_byt=[zeros(lenH-1,1);lenTyp];  %all Peri
end

vit=[];
for j=1:length(H)
    vit=[vit; j*ones(n_byt(j),1)];
end

gVal=input('g value 2=.5nS, 22=1.25nS, 3=2nS, 33=3nS, 4=4nS)? : ');
sVal=input('Gradient in G (1:P=15xC, 2:C=15XP, 3=same,.5gap])? : ');
switch gVal
    case 2
        gvf=0.0005;
    case 22
        gvf=0.00125;
    case 3
        gvf=0.002;
    case 33
        gvf=0.003;
    case 4
        gvf=0.004;
end
switch sVal
    case 1
        [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,0,ind_j); %MAKE cells near peripheral STRONGER gap
        G_cup=G_cup*gvf;
    case 2
        [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j); %MAKE cells near CENTER STRONGER gap
        G_cup=G_cup*gvf;
    case 3
        [G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j); 
        G_cup(G_cup~=0)=0.5; %MAKE cells have same g_param
        G_cup=G_cup*gvf;
end
%calc effective g for each of the 9 oscill; use Het 'type' even for Homo!
%  because using same ind_j from H_dataHexGrd.mat
ampB=zeros(lenTyp,1);
mapType=[(1:2:15)'; 18]; %assuming lenType=9
tot_G=sum(G_cup,2);
for j=1:lenTyp
    ampB(j,1)=mean(tot_G(ind_j==mapType(j)));
end

load ../disi

ind_om=1./isi_L(vit); %get individual frequencies
ind_om=ind_om./ind_om(1); %scale so slowest (inner) is 1

load('../Parms_c2p','C') %get conductances

phs_v=cell(lenH,1);
Hsc=cell(lenH,1);
ampA=zeros(lenH,1); %need all amp, for scaling below
ampA(mapType)=ampB; %only non-zero if used

for j=1:lenH
    %ALL of them (heterogen)
    phs_v{j,1}=tmh{j}./tmh{j}(end); %scaled so between [0,1)
    Hsc{j,1}=H{j}./C(j)*ampA(j);     %scaled H-fcn (didn't use C in XPP)
end


Lt=10000;
dt=0.001;
t=(dt:dt:dt*Lt)';

if(exist('thet'))
    thet(:,1)=thet(:,end);
else
    thet(:,1)=rand(lenTyp,1);
end


for j=2:Lt
    th_prev=thet(:,j-1);
    
    thet(1,j)=th_prev(1)+dt*(ind_om(1)+...
        interp1q(phs_v{vit(1)},Hsc{vit(1)},mod(th_prev(2)-th_prev(1),1)));
    for k=2:(lenTyp-1)
        thet(k,j)=th_prev(k)+dt*(ind_om(k)+...
            +interp1q(phs_v{vit(k)},Hsc{vit(k)},mod(th_prev(k+1)-th_prev(k),1))...
            +interp1q(phs_v{vit(k)},Hsc{vit(k)},mod(th_prev(k-1)-th_prev(k),1)) );
    end
    thet(lenTyp,j)=th_prev(lenTyp)+dt*(ind_om(end)+...
        interp1q(phs_v{vit(lenTyp)},Hsc{vit(lenTyp)},mod(th_prev(lenTyp-1)-th_prev(lenTyp),1)));

end

OP=abs( sum(exp(2*pi*sqrt(-1)*thet))./lenTyp ); %order param


xd=diff(mod(thet',1))';
ySpk=1*(xd<0);
figure
spy(ySpk,24)
axis square

figure
hold on
caxis([0 1]);
pcolor(repmat(t',lenTyp,1),repmat((1:lenTyp)',1,Lt),mod(thet,1))
shading interp
colorbar
set(gca,'FontSize',18)
xlabel('Time')
ylabel('Neuron Index')

figure
plot(t,OP)
set(gca,'FontSize',18)
ylabel('Order Parameter')
xlabel('Time')