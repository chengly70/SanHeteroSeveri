
%Script to plot characteristic time-scale (with theory)

lenH=18; %assuming lenN=lenH
lenTyp=9;
mapType=[(1:2:15)'; 18]; %assuming lenType=9

load ../disi %for scaling time

load svGavgCompHex %getting average g (gAvg_comp) from full sims (calc_gAvg_fromFullHex.m)

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
CharTm_ss=cell(3,3); 
LrgNzEig=cell(3,3);
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
                scl_time=isi_L(1);
                
            case 2
                n_byt=[lenTyp;zeros(lenH-1,1)];  %all Cent
                scl_time=isi_L(1);
                
            case 3
                n_byt=[zeros(lenH-1,1);lenTyp];  %all Peri
                scl_time=isi_L(lenH);
                
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
            
            Omga=prs(1);
            oxi=[0;prs(2:end)];

            load Hfs %the uncoupled H's, must be the same size as isi_L
            
            load('../Parms_c2p','C') %get conductances

            phs_v=cell(lenTyp,1);
            Hsc=cell(lenTyp,1);

            for j=1:lenTyp
                %not all Hsc, only needed; 
                phs_v{j,1}=tmh{vit(j)}./tmh{vit(j)}(end); %scaled so between [0,1)
                Hsc{j,1}=H{vit(j)}./C(vit(j))*ampB(j);     %scaled H-fcn (ampB from loaded sHex_[].mat)
            end

            %--- fund. matrix ---
            A_p=zeros(lenTyp,lenTyp);
            h_drv=0.01; %step for H derivative
            
            % get avg H fcns on corase grid
            phCrs=(0:0.01:1)';
            Havg=zeros(length(phCrs),lenTyp);
            for j=1:lenTyp
                for k=1:length(phCrs)
                    Havg(k,j)=interp1q(phs_v{j},Hsc{j},phCrs(k));
                end
            end
            rgh_p=interp1q(phCrs,Havg(:,1),mod(oxi(2)-oxi(1)+h_drv,1));
            lft_p=interp1q(phCrs,Havg(:,1),mod(oxi(2)-oxi(1)-h_drv,1));
            derv_p=(rgh_p-lft_p)/h_drv;
            A_p(1,1:2)=A_p(1,1:2)+[-derv_p derv_p];
            for j=2:(lenTyp-1)
                rgh_p=interp1q(phCrs,Havg(:,j-1),mod(oxi(j-1)-oxi(j)+h_drv,1));
                lft_p=interp1q(phCrs,Havg(:,j-1),mod(oxi(j-1)-oxi(j)-h_drv,1));
                derv_p=(rgh_p-lft_p)/h_drv;
                A_p(j,j-1)=A_p(j,j-1)+derv_p;
                A_p(j,j)=A_p(j,j)-derv_p;
                
                rgh_p=interp1q(phCrs,Havg(:,j),mod(oxi(j+1)-oxi(j)+h_drv,1));
                lft_p=interp1q(phCrs,Havg(:,j),mod(oxi(j+1)-oxi(j)-h_drv,1));
                derv_p=(rgh_p-lft_p)/h_drv;
                A_p(j,j)=A_p(j,j)-derv_p;
                A_p(j,j+1)=A_p(j,j+1)+derv_p;
            end
            rgh_p=interp1q(phCrs,Havg(:,lenTyp),mod(oxi(lenTyp-1)-oxi(lenTyp)+h_drv,1));
            lft_p=interp1q(phCrs,Havg(:,lenTyp),mod(oxi(lenTyp-1)-oxi(lenTyp)-h_drv,1));
            derv_p=(rgh_p-lft_p)/h_drv;
            A_p(lenTyp,lenTyp-1:lenTyp)=A_p(lenTyp,lenTyp-1:lenTyp)+[derv_p -derv_p];

            [V,D]=eig(A_p);
            [eig_vals,indEig]=sort(diag(D),'descend');
            
            LrgNzEig{svl_v(sInd),hethom_v(hInd)}(gInd,1)=eig_vals(2);
            CharTm_ss{svl_v(sInd),hethom_v(hInd)}(gInd,1)=scl_time/abs(eig_vals(2));
            
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
        plot(gVal_nS{svl_v(sInd),hethom_v(hInd)},CharTm_ss{svl_v(sInd),hethom_v(hInd)},...
            '-','LineWidth',2,'color',cc(hInd,:))
    end
    legend('Het','Hom (C)','Hom (P)')
    set(gca,'FontSize',18)
    box off
    xlabel('Gap Junction Strength (nS)')
    ylabel('Characteristic Time-scale (s)')
end


