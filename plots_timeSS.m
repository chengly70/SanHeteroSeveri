
%script to calc the time to reach SS (traveling wave) FULL sims
% calc range of voltage for each cell to calc Tolerance
%   epsTlf = epsTlf*rangVlt*(# pts); epsSynch

%using previously calc Freq of trav wave to see when SS reached
load Freq_fullModel 

load('dhex_het_io_g22_s3.mat','Nz') %load any het file to get Nz

epsTlf=0.001; %used epsTlf=0.0017 for het ONLY

rangVlt=zeros(Nz,1); %update voltageRange for each sim

load thry/svGavgCompHex %getting average g (gAvg_comp) from full sims (calc_gAvg_fromFullHex.m)

%vectors to denote values to range over
%svl_v=[1;2;3];    %1:P=15XC, 2:C=15XP, 3=same,.5*g_gap
%hethom_v=[1;2;3]; %1=het, 2=homC, 3=homP
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

cc=[zeros(1,3); 0 0 1; 0 .5 0]; %color-code for het/hom: Het=Black, homC=blue, homP=green

%variables to be saved
Time_ss=cell(3,3); 
%first col is valid g for heterog; 2nd col is all g's (hom C & P)
%first row is S1, second row is S2, 3rd row is S3

t1=(0:0.001:20)'; %finer time mesh
dt=t1(2)-t1(1);
% --- starting the large loop to loop through everything ---
for sInd=1:3
    
    for hInd=1:3
        
        if(hInd==1) %special cases
            epsTlf=0.0017;
        else
            epsTlf=0.001; %set to orig value (line 12 above)
        end
        
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
            
            % -- calc of timeSS from large sims ---
            load(flname)
            
            %only saving voltages so don't have to set (:,1:Nz)
            volt_pop=interp1(time,V_out,t1,'pchip'); %interp on finer time mesh
            
            rangVlt=( max(volt_pop) - min(volt_pop) )'; %range of voltage

            freqTw=Frq_SA_ss{svl_v(sInd),hethom_v(hInd)}(gInd); %prev calc Freq
            szVecMt=round(1/freqTw/dt);
            
            Err=100*ones(Nz,1); %check on cell-cell basis
            j=szVecMt+1;
            while(sum(Err>epsTlf*szVecMt*rangVlt)>0 && j<length(volt_pop)-2-szVecMt)
                v_to_match=volt_pop(j+1:j+1+szVecMt,:); %update to next cycle
                v_to_match1=volt_pop(j:j+szVecMt,:); 
                v_to_match2=volt_pop(j+2:j+2+szVecMt,:); 
                for k=1:Nz %use vecnorm for 2017b & higher
                    Err(k,1)=min([norm(volt_pop(j-szVecMt:j,k)-v_to_match(:,k)) ...
                    norm(volt_pop(j-szVecMt:j,k)-v_to_match1(:,k)) ...
                        norm(volt_pop(j-szVecMt:j,k)-v_to_match2(:,k))]);
                end
                j=j+1;
            end
            %save calc time here
            Time_ss{svl_v(sInd),hethom_v(hInd)}(gInd,1)=t1(j-szVecMt);

        end
        
        %plot here;
        figure(sInd)
        plot(gVal_nS{svl_v(sInd),hethom_v(hInd)},Time_ss{svl_v(sInd),hethom_v(hInd)},...
            '*-','MarkerSize',18,'color',cc(hInd,:))
    end
    set(gca,'FontSize',20)
    box off
    
end

save dTimeSS_all gVal_nS svl_v hethom_v Time_ss epsTlf