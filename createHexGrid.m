N = 10;  % should be even, but can be adjusted for odd number
[X,Y] = meshgrid(-(N-1):N-1);
X(1:2:end,:) = X(1:2:end,:) + .5;

dx = 0.5;
delta = dx;

if mod(N,2)  % if odd number of layers
    X(N,end) = nan; 
    Y(N,end) = nan;
end
    
for i = 1:N-1
    xtmp = X(N-i,:);
    k = (xtmp<-(N-1)+delta) | (xtmp>(N-1)-delta);
    X(N-i,k) = nan;
    Y(N-i,k) = nan;
    
    xtmp = X(N+i,:);
    k = (xtmp<-(N-1)+delta) | (xtmp>(N-1)-delta);
    X(N+i,k) = nan;
    Y(N+i,k) = nan;
    delta = delta + dx;
end
Y = Y*sqrt(3)/2;
 x = X(:); y = Y(:);
 i = ~isnan(x);
x = x(i); y = y(i);
d = sqrt(x.^2 + y.^2);

% figure
% plot(x,y,'o');

%so variables match
Nz=length(x);
X_l=x;
Y_l=y;

%vector to identify different cell types
ind_j=zeros(Nz,1);
vc_Nz=(1:Nz)';
%set cell at origin to 'Center'; also 6 nearest neighbors
[m,idC]=min(d); 
ind_j(idC)=1;
idC = (d>0.99).*(d<1.01).*(ind_j==0);  %cells within distance 1 & NOT assigned
ind_j(idC==1)=1;
for j=2:(N-1)
    preV=vc_Nz(ind_j==(j-1));
    for iPre=1:length(preV)
        d_to_Pre=sqrt((X_l-X_l(preV(iPre))).^2+(Y_l-Y_l(preV(iPre))).^2);
        idC = (d_to_Pre>0.99).*(d_to_Pre<1.01).*(ind_j==0);
   
        ind_j(idC==1)=j;
    end
end

%set 2->3, 3->5, 2*ind_j-1,.. so spans 1:2:17
for j=(N-1):-1:2
    ind_j(ind_j==j) = 2*j-1;
end
%set #17 to #18 to compare to homP
ind_j(ind_j==17)=18;

%save some variables
save H_dataHexGrd Nz X_l Y_l ind_j

