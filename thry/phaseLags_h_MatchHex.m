%numerically calculate the ensemble phase Omga and offsets using fminunc
% for HETEROG (C<->P) pacemaker cells
% similar to phaseLags_het.m BUT allowing different distrib (#'s of cell
% types)
% !!ASSUMING we have: lenH, lenTyp, ampB, y0 from sim_phaseOscill_Match.m & get_y0_fromSims.m
%  as well as hetOrHom, gVal, sVal, ind_om, vit, n_byt

load Hfs %the uncoupled H's, must be the same size as isi_L

load ../disi %get isi_L

load('../Parms_c2p','C') %get conductances

phs_v=cell(lenH,1);
Hsc=cell(lenH,1);

for j=1:lenH 
    %all heterogen
    phs_v{j,1}=tmh{j}./tmh{j}(end); %scaled so between [0,1)
    Hsc{j,1}=H{j}./C(j)*ampA(j);     %scaled H-fcn
end

%start iteration to solve transcendental eqns self-consistently

%initial conditions

prs0 = y0; %from running get_y0_fromSims

fxn = @(prs)eqn_distr(prs,ind_om,Hsc,phs_v,vit);
A = [];
b = [];
Aeq = [];
beq = [];
lb=zeros(lenTyp,1);                  %lower bound
ub=[max(ind_om) ; ones(lenTyp-1,1)]; %upper bound


opts = optimset('MaxFunEvals', 500000, 'MaxIter', 50000,'Algorithm','sqp');%, 'Display', 'off');
%exitflag = 0;
%maxNumRuns = 1;
[prs,fval,exitflag,oup,lams,grd,hess] = fmincon(fxn,prs0,A,b,Aeq,beq,lb,ub,[],opts);

fval

[obj,ob_pcs]=eqn_distr(prs,ind_om,Hsc,phs_v,vit);

save(['sHex_H',num2str(hetOrHom),'_g',num2str(gVal),'_s',num2str(sVal)],'prs','prs0','fval','ob_pcs','ampB','y0','lenTyp','vit')

Omga=prs(1);
oxi=[0;prs(2:end)];

t=(0:.05:20)';
Lt=length(t);
X_app=zeros(lenTyp,Lt); %approx to phase model
X_app=mod(Omga*repmat(t',lenTyp,1)+repmat(oxi,1,Lt),1);


figure
hold on
caxis([0 1]);
pcolor(repmat(t',lenTyp,1),repmat((1:lenTyp)',1,Lt),X_app)
shading interp
colorbar
set(gca,'FontSize',18)
xlabel('Time')
ylabel('Neuron Index')

disp_on=1;

if(disp_on)
    
% show deriv of Hj at soln
H_atX=zeros(lenTyp,2);
% get avg H fcns on corase grid
phCrs=(0:0.01:1)';
Havg=zeros(length(phCrs),lenTyp);
for j=1:lenTyp
    for k=1:length(phCrs)
    Havg(k,j)=interp1q(phs_v{vit(j)},Hsc{vit(j)},phCrs(k));
    end
end
figure
hold on
cc=copper(lenTyp);
for j=1:lenTyp
    plot(phCrs,Havg(:,j),'color',cc(j,:),'LineWidth',2)
end
for j=1:lenTyp
    if(j>1)

    H_atX(j,1)=interp1q(phs_v{vit(j)},Hsc{vit(j)},mod(oxi(j-1)-oxi(j),1));
        plot(mod(oxi(j-1)-oxi(j),1),H_atX(j,1),'b*','MarkerSize',14)
    end
    if(j<lenTyp)

    H_atX(j,2)=interp1q(phs_v{vit(j)},Hsc{vit(j)},mod(oxi(j+1)-oxi(j),1));
        plot(mod(oxi(j+1)-oxi(j),1),H_atX(j,2),'r*','MarkerSize',14)
    end
end
set(gca,'FontSize',18)
xlabel('Phase')

end
