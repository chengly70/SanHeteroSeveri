
%Script to plot frequencies from reduced theory

lenH=18; %assuming lenN=lenH
lenTyp=9;

load ../disi %for scaling time

load svGavgCompHex %getting average g (gAvg_comp) from full sims (calc_gAvg_fromFull.m)

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

%variables to be saved
gVal_nS=cell(3,3); %actual values in full sims(nano-Siemens), for plotting
FrqNetw_ss=cell(3,3); 
%first col is valid g for heterog; 2nd col is all g's (hom C & P)
%first row is S1, second row is S2, 3rd row is S3

cc=[zeros(1,3); 0 0 1; 0 .5 0]; %color-code for het/hom: Het=Black, homC=blue, homP=green

% --- starting the large loop to loop through everything ---
for sInd=1:3

    figure
    hold on
    
    for hInd=1:3
        switch hInd
            case 1
                n_byt=[repmat([1;0],8,1); 0; 1]; %het type: 1,3,5,..,15,18
                scl_Frq=1/isi_L(1);
            case 2
                n_byt=[lenTyp;zeros(lenH-1,1)];  %all Cent
                scl_Frq=1/isi_L(1);
            case 3
                n_byt=[zeros(lenH-1,1);lenTyp];  %all Peri
                scl_Frq=1/isi_L(lenH);
        end
        vit=[];
        for j=1:lenH
            vit=[vit; j*ones(n_byt(j),1)];
        end
            
        for gInd=1:length(gval_c{sInd,hInd})
            flname='sHex_H';
            flname=[flname,num2str(hethom_v(hInd)),'_g'];
            
            gVal=gval_c{svl_v(sInd),hethom_v(hInd)}(gInd);
            flname=[flname,num2str(gVal),'_s'];
            
            flname=[flname,num2str(svl_v(sInd))];
            
            load(flname)
            
            %only use largest one; checked that the others have smallER dot
            %products; 
            FrqNetw_ss{svl_v(sInd),hethom_v(hInd)}(gInd,1)=prs(1);
            
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
        % scl_Time is to scale by CL of 1 or 18 (P), see above
        plot(gVal_nS{svl_v(sInd),hethom_v(hInd)},scl_Frq*FrqNetw_ss{svl_v(sInd),hethom_v(hInd)},...
            '-','LineWidth',2,'color',cc(hInd,:))
    end
    legend('Het','Hom (C)','Hom (P)')
    set(gca,'FontSize',20)
    box off
    xlabel('Gap Junction Strength (nS)')
    ylabel('Wave Frequency (1/s)')
end

