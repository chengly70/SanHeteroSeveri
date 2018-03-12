%script to get simulated 'periodic' steady-state 
% AFTER running sim_phaseOscill_Match.m; assuming have: t, dt, thet

if(exist('lenN'))
    lny=lenN;
else
    lny=9; %for Hex grid
end

y0=zeros(lny,1); %output

tm_start=input('Enter time when steady-state is reached: ');

id_Strt=round( (tm_start-t(1))/dt + 1);

stVar=mod(thet(:,id_Strt:end),1);

numInst=size(stVar,2); %number of instances

ofx=zeros(lny-1,1); %offsets, xi_j, relative to theta(1)

ofx=mean( mod(stVar(2:end,:)-repmat(stVar(1,:),lny-1,1),1), 2);
y0(2:end)=ofx;


%power spectrum calc (and auto corr)
mLag=3; %in sec
mLag_i=floor(mLag/dt);
t_aut=(-mLag:dt:mLag)';
ind_z=round((0-t(1))./dt+1); %index of t=0

Fs = 1/dt; % Sampling frequency
T = dt; % Sample time

for j=1:lny
    Acr(:,j)=xcorr(stVar(j,:),mLag_i,'unbiased'); %autcorr for stim (raw)
    Acr(:,j)=Acr(:,j)-mean(stVar(j,:))^2;
    
    L = length(Acr(:,j)); % Length of signal
    xdft = fft(Acr(:,j));
    Pss(:,j) = 1/(L*Fs)*abs(xdft(1:floor(L/2)+1)).^2;
end
%get max freq for each cell:
[mx,f_id]=max(Pss);
freq = 0:Fs/L:Fs/2;

%set y0(1) to mean of max freq:
y0(1)=mean(freq(f_id));

if(0)
    cc=jet(lny);
    figure(13)
    hold on
    figure(14)
    hold on
    for j=1:lny
        figure(13) %autocorre
        plot(t_aut,Acr(:,j),'color',cc(j,:),'LineWidth',1.5)
        
        figure(14)
        plot(freq,Pss(:,j),'color',cc(j,:),'LineWidth',1.5)
    end
end