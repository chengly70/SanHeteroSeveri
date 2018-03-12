
%plots spatial location AND gray segments to denote connectivity
%with Hexagonal Grid.  assume file H_dataHexGrd exists
%!! to color-code cell types, MUST have Bind=Nz x 1 variable that has cell
%type (1,2,..,18) in each entry 
%Can get Bind from run_[].m scripts, see code

load H_dataHexGrd %(X_l,Y_l), ind_j, Nz

gp_lim=1.05; %the Euclidean length limit for which there is no gap junction
            %coupling; units are arbitrary X_l,Y_l, see ax/ay in ellpse_coor.m
[G_cup]=nn_coupVryStrHex(X_l,Y_l,gp_lim,1,ind_j); %MAKE cells near CENTER STRONGER gap


%******************Initial Conditions********************, In->Out T.W.
sclDist=sqrt(X_l.^2+Y_l.^2);
sclDist=sclDist./max(sclDist);
col_IC=1-sclDist; %color-codes Init Cond (how close/far to spiking)

figure
hold on
axis square

Ns=length(X_l);
for j=1:Ns
    for k=j+1:Ns
        if(G_cup(j,k)~=0)
            plot(X_l([j;k]),Y_l([j;k]),'color',.5*ones(1,3))
        end
    end
end
    

colormap('jet')
scatter(X_l,Y_l,50,col_IC,'fill') % ind_j identifies cell type
colorbar
axis off