
%script to show Power Spects of Volt for dhex_het_io_g2_s1.mat, as well 
%as the Voltage Traject (and location on 2D hex, see commented out snippet)

load thry/svGavgCompHex %getting average g (gAvg_comp) from full sims (calc_gAvg_fromFull.m)

minFreq=2; %min freq in Hz that traveling wave is allowed to be

%variables to be saved
gVal_nS=cell(3,3); %actual values in full sims(nano-Siemens), for plotting
Frq_SA_ss=cell(3,3); 
%first col is valid g for heterog; 2nd col is all g's (hom C & P)
%first row is S1, second row is S2, 3rd row is S3

t1=(0:0.001:20)'; %finer time mesh

% Specficy dhex_het_io_g2_s1.mat with vars: sInd, hInd, gVal
sInd=1;
hInd=1;
gVal=2; %from [2,22,3,33,4]

switch hInd
    case 1
        flname='dhex_het_';
    case 2
        flname='dhexhomC_';
    case 3
        flname='dhexhomP_';
end
flname=[flname,'io_']; %assuming all in->out IC
flname=[flname,'g',num2str(gVal),'_s',num2str(sInd)];

% -- calc of freq from large sims ---
load(flname)
V_outF=interp1(time,V_out,t1,'pchip'); %interp on finer time mesh
volt_pop=V_outF;%only saving voltages so don't have to set (:,1:Nz)
v_mx=max(max(volt_pop));
v_mn=min(min(volt_pop));
%for calculating effective phase
dt=t1(2)-t1(1);
tm_start=10;
id_str=round(tm_start/dt+1); %approx index for when time=2sec; assume transient settled after
v_mxO=(max(volt_pop(id_str:end,:)))';
v_mnO=(min(volt_pop(id_str:end,:)))';
len_op=length(t1)-id_str+1;
phas_var=(volt_pop(id_str:end,:)-repmat(v_mnO',len_op,1))./(repmat(v_mxO',len_op,1)-repmat(v_mnO',len_op,1));
OrdP=abs( sum(exp(2*pi*sqrt(-1)*phas_var'))./Nz )';
tmO=tm_start+(0: (t1(end)-tm_start)/(len_op-1) : (t1(end)-tm_start) )';

mLag=10; %in sec
mLag_i=floor(mLag/dt);
t_aut=(-mLag:dt:mLag)';
Fs = 1/dt; % Sampling frequency
T = dt; % Sample time
for j=1:Nz
    Acr(:,j)=xcorr(phas_var(:,j),mLag_i,'unbiased'); %autcorr for stim (raw)
    Acr(:,j)=Acr(:,j)-mean(phas_var(:,j))^2;
    
    L = length(Acr(:,j)); % Length of signal
    xdft = fft(Acr(:,j));
    Pss(:,j) = 1/(L*Fs)*abs(xdft(1:floor(L/2)+1)).^2;
end
%get max freq for each cell:
freq = 0:Fs/L:Fs/2;
stInd=round(minFreq*L/Fs)+1; %only freq above 2
[mx,f_id]=max(Pss(stInd:end,:));
f_id=f_id+stInd-1; %adjust b/c min freq offset


%plot the power spectrums, show it is multimodal    
thres_ps=.004; %threshold for determining multimodal (!!!DEPENDS on params!!!)
cells_mm=[]; %index of cells that have multimodal Power Spects
cel_sing=[]; %index of cells that have a single peak Power Spects
figure
hold on
for j=1:Nz
    if(sum(Pss(:,j)>thres_ps)>1 && j==141) %multi-modal
        cells_mm=[cells_mm; j]; 
        plot(freq,Pss(:,j),'m','LineWidth',2)
    elseif(sum(Pss(:,j)>thres_ps)>1)
        cells_mm=[cells_mm; j]; 
        plot(freq,Pss(:,j),'b','LineWidth',2)
    else
        plot(freq,Pss(:,j),'k','LineWidth',1)
        cel_sing=[cel_sing; j];
    end
end
set(gca,'FontSize',20)
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (mV^2)')
set(gca,'XLim',[2 6])

%plot voltage
dist_C=sqrt(X_l.^2+Y_l.^2); %distance to center (0,0)
[dval,ind_inOut]=sort(dist_C);
cc=copper(Nz);
figure
hold on
for j=1:length(cel_sing) %want these behind
    if(ismember(ind_inOut(j),cel_sing))
        plot(t1,volt_pop(:,ind_inOut(j)),'color',cc(j,:))
    end
end
for j=1:length(cells_mm)
    if(ismember(ind_inOut(j),cells_mm) && ind_inOut(j)==141)
        plot(t1,volt_pop(:,ind_inOut(j)),'m','LineWidth',2)
    elseif(ismember(ind_inOut(j),cells_mm))
        plot(t1,volt_pop(:,ind_inOut(j)),'b','LineWidth',2)
    end
end
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Voltage (mV)')
set(gca,'YLim',[-70 50])


%Panel C)
open 2Dmap_Hex.fig    
for j=1:length(cells_mm)
    plot(X_l(cells_mm(j)),Y_l(cells_mm(j)),'b.','MarkerSize',30)
end
plot(X_l(141),Y_l(141),'c.','MarkerSize',30) %this is the weird one
