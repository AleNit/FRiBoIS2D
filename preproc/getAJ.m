
% compute area and polar moment of are w.r.t. the rotation center

function [AA,JJ]=getAJ(z,y,c)

% distribute points inside the curve
z1=linspace(min(z)-0.1,max(z)+0.1,60);
y1=linspace(min(y)-0.1,max(y)+0.1,30);
[Z1,Y1]=meshgrid(z1,y1);
Z1=reshape(Z1,[],1);
Y1=reshape(Y1,[],1);
in = inpolygon(Z1,Y1,z,y);
zp=Z1(in);
yp=Y1(in);

zp=[zp;z];
yp=[yp;y];

% get triangulation
P=[zp,yp];
DT = delaunayTriangulation(P);

tc=incenter(DT);

% compute triangle area
np=length(DT.ConnectivityList(:,1));
A=zeros(np,1);
for i=1:np
     x1 = DT.Points(DT.ConnectivityList(i,:),:);
     A(i) = 1/2.*abs((x1(2,1)-x1(1,1)).*(x1(3,2)-x1(1,2))-(x1(3,1)-x1(1,1)).*(x1(2,2)-x1(1,2)));     
end

AA=sum(A);
JJ=sum( ( (tc(:,1)-c(1)).^2 + (tc(:,2)-c(2)).^2 ).*A );

% debug plot
% figure(10)
% triplot(DT)
% hold on
% plot(tc(:,1),tc(:,2),'xr')
% axis equal
% xlabel('z/D','FontSize',14)
% ylabel('y/D','FontSize',14)
% legend('triangulation','centers','FontSize',14)


end