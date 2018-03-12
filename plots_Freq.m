
%script to calc the freq of traveling waves for FULL sims
%assuming we do not have steady state until 10s (conservative!)

load thry/svGavgCompHex %getting average g (gAvg_comp) from full sims (calc_gAvg_fromFull.m)

%vectors to denote values to range over
svl_v=[1;2;3];    %1:P=15XC, 2:C=15XP, 3=same,.5*g_gap
hethom_v=[1;2;3]; %1=het, 2=homC, 3=homP
%first col is valid g for heterog; 2nd col is all g's (hom C & P)
gval_c=cell(3,3); %first row is S1, second row is S2, 3rd row is S3
gval_c{1,1}=[22;3;33;4]; %found & checked these manually
gval_c{2,1}=[2;22;3;33;4];
gval_c{3,1}=[2;22;3;33;4];
for j=1:3
    for k=2:3
        gval_c{j,k}=[2;22;3;33;4]; %all are ok..?
    end
end

minFreq=2; %min freq in Hz that traveling wave is allowed to be

cc=[zeros(1,3); 0 0 1; 0 .5 0]; %color-code for het/hom: Het=Black, homC=blue, homP=green

%variables to be saved
gVal_nS=cell(3,3); %actual values in full sims(nano-Siemens), for plotting
Frq_SA_ss=cell(3,3); 
%first col is valid g for heterog; 2nd col is all g's (hom C & P)
%first row is S1, second row is S2, 3rd row is S3

t1=(0:0.001:20)'; %finer time mesh
% --- starting the large loop to loop through everything ---
for sInd=1:3
    figure(sInd) %assuming Figs1,2,3 already open
    for hInd=1:3
        
        for gInd=1:length(gval_c{sInd,hInd})
            
            switch hInd
                case 1
                    flname='dhex_het_';
                case 2
                    flname='dhexhomC_';
                case 3
                    flname='dhexhomP_';
            end
            flname=[flname,'io_']; %assuming all in->out IC
            gVal=gval_c{svl_v(sInd),hethom_v(hInd)}(gInd);
            flname=[flname,'g',num2str(gVal),'_s',num2str(svl_v(sInd))];
            
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

            %save calc freq here
            Frq_SA_ss{svl_v(sInd),hethom_v(hInd)}(gInd,1)=mean(freq(f_id));

            %get gAvg value for x-axis to plot
            switch gVal
                case 2
                    gVal_nS{svl_v(sInd),hethom_v(hInd)}(gInd,1)=gAvg_comp(sInd,1);
                case 22
                    gVal_nS{svl_v(sInd),hethom_v(hInd)}(gInd,1)=gAvg_comp(sInd,2);
                case 3
                    gVal_nS{svl_v(sInd),hethom_v(hInd)}(gInd,1)=gAvg_comp(sInd,3);
                case 33
                    gVal_nS{svl_v(sInd),hethom_v(hInd)}(gInd,1)=gAvg_comp(sInd,4);
                case 4
                    gVal_nS{svl_v(sInd),hethom_v(hInd)}(gInd,1)=gAvg_comp(sInd,5);
            end
        end
        
        %plot here;
        plot(gVal_nS{svl_v(sInd),hethom_v(hInd)},Frq_SA_ss{svl_v(sInd),hethom_v(hInd)},...
            '*','MarkerSize',18,'color',cc(hInd,:))
    end
end

save Freq_fullModel gVal_nS svl_v hethom_v Frq_SA_ss