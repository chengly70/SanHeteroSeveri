%script to plot voltage traces only; time goes to 20s but zoom-in to
%actually see it

hetOrHom=input('Heterog (1) or Homog (2=Cent, 3=Periph)? : ');
switch hetOrHom
    case 1
        flname='dhex_het_';
    case 2
        flname='dhexhomC_';
    case 3
        flname='dhexhomP_';
end
flname=[flname,'io_'];

gVal=input('g value (2=.5nS, 22=1.25nS, 3=2nS, 33=3nS, 4=4nS)? : ');
sVal=input('Gradient in G (1:P=15xC, 2:C=15XP, 3=same,.5gap])? : ');
flname=[flname,'g',num2str(gVal),'_s',num2str(sVal)];

load(flname) %get large-scale sim

volt_pop=V_out;%only saving voltages so don't have to set (:,1:Nz)

v_mx=max(max(volt_pop));
v_mn=min(min(volt_pop));

%show all voltage traj
cc=copper(Nz);
figure
hold on
dist_C=sqrt(X_l.^2+Y_l.^2); %distance to center (0,0)
[dval,ind_inOut]=sort(dist_C);

for j=1:Nz
    plot(time,volt_pop(:,ind_inOut(j)),'color',cc(j,:))
end
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Voltage (mV)')
if(hetOrHom==2) %center cells only
    axis([0 5 -60 30])
else
    axis([0 5 -70 50]) %has at least some peripheral cells
end
