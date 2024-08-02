
% create input files for the galloping rectangle test
% Choose here grid setting and build geometry. The input files
% are taken from the ./input_template folder and modified accordingly

% Alessandro Nitti, Polytechnic University of Bari (2024)

clc 
clear
close all
kf=1;


%% input parameters
cname='test_GallRect1DOF';      % case name
wrt=false(1);                   % create input folder
l2e=60.0;                       % domain extension in the y direction
l3e=40.0;                       % domain extension in the z direction
nD=50;                          % number of Eulerian nodes over the short rectangle side
z0=20;                          % z coordinate of the body rotation center
y0=ceil(l2e/2);                 % y coordinate of the body rotation center
rr=0.5;                         % lagrangian to eulerian grid spacing ratio 



%% generate cartesian grid
int3=[round(z0-3),round(z0+7)];
nint3=[ceil(nD/2),(int3(2)-int3(1))*nD,ceil(nD*4)+6];
x3n=cubicstr(int3,nint3,0,l3e);
n3=sum(nint3);

int2=[ceil(l2e/2)-5,ceil(l2e/2)+5];
nint2=[ceil(0.6*nD),(int2(2)-int2(1))*nD,ceil(0.6*nD)+1];        
x2n=cubicstr(int2,nint2,0,l2e); 
n2=sum(nint2);

dymin=10^6;
for i=2:n2    
    tmp=x2n(i)-x2n(i-1);
    dymin=min(dymin,tmp);
end
dzmin=10^6;
for i=2:n3    
    tmp=x3n(i)-x3n(i-1);
    dzmin=min(dzmin,tmp);
end
dmin=min(dzmin,dymin);

disp(['n2 = ' num2str(n2)])
disp(['n3 = ' num2str(n3)])
disp(['total number of points: ' num2str(n2*n3)])

[X3n,Y3n]=meshgrid(x3n,x2n);
Z3n=zeros(n2,n3);

figure(kf); kf=kf+1;
surf(X3n,Y3n,Z3n,'facecolor','w')
hold on
view(2)
axis equal
axis([0,l3e,0,l2e])


%% create geometry
L=3;
lp=-L/2:mean([dzmin,dymin])*rr:L/2;
sp=-0.5:mean([dzmin,dymin])*rr:0.5;
zcp=[z0+lp, (z0+L/2).*ones(1,length(sp)-2), flip(z0+lp), (z0-L/2).*ones(1,length(sp)-2)]';
ycp=[(y0-0.5).*ones(1,length(lp)), sp(2:end-1)+y0, (y0+0.5).*ones(1,length(lp)), flip(sp(2:end-1))+y0]';
nlm=length(zcp)-1;

% compute area and polar moment of are w.r.t. the rotation center
[AA,JJ]=getAJ(zcp,ycp,[z0,y0]);
disp(['area: ',num2str(AA)])
disp(['polar moment of area: ',num2str(JJ)])


%% create geometrical features for IB treatment
[zlm,ylm,Alm,nor]=getlmn(nlm,zcp,ycp,-1,dzmin,dymin);

scatter(zcp,ycp,1,'markeredgecolor','b')
scatter(zlm,ylm,3,'markeredgecolor','r')
scatter(zlm(1),ylm(1),10,'filled')
scatter(zlm(4),ylm(4),10,'filled')

nlm=length(zlm);
for i=1:nlm
    yplt=[ylm(i),nor(i,1)];
    zplt=[zlm(i),nor(i,2)];    
    xplt=[0.0,0.0];
    plot3(zplt,yplt,xplt ,'color','r')
    hold on
end

xlabel('z');ylabel('y')

disp(['number of Lagrangian markers: ' num2str(nlm)])


%% Write data for simulation input
if (wrt)

    % create input folder
    infold=strcat('../',cname,'/input_FSI/');
    if (~exist(infold, 'dir'))
        mkdir(infold)
    end


    % grid ccordinates    
    name=strcat(infold,'zcoord.dat');
    fID = fopen(name,'w');
    for i=1:n3-1
        fprintf(fID,'%16.14f \n',x3n(i));    
    end
    fprintf(fID,'%16.14f',x3n(end));    
    fclose(fID);   

    name=strcat(infold,'ycoord.dat');
    fID = fopen(name,'w');
    for i=1:n2-1   
        fprintf(fID,'%16.14f \n',x2n(i));    
    end
    fprintf(fID,'%16.14f',x2n(end));    
    fclose(fID);   
    

    % edge vertices 
    name=strcat(infold,'lmdata.in');
    fID = fopen(name,'w');
    for i=1:nlm
        prt=[ycp(i),zcp(i)];
        fprintf(fID,'%16.14f \t',prt);
        fprintf(fID,'\n');
    end
    fprintf(fID,'%16.14f \t',[ycp(end),zcp(end)]);
    fclose(fID);


    % create input files form template: rigid_par.in
    filen='./input_template/rigid_par.in';
    fileID=fopen(filen);        
    FC=textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);

    str=sprintf('%05.4f, %05.4f',y0,z0);
    FC{1}{3}(1:length(str))=str;

    fID=fopen(strcat(infold,'rigid_par.in'),'w');
    for ii=1:length(FC{1})
        fprintf(fID,'%s \n',FC{1}{ii});
    end
    fclose(fID);


    % create input files form template: fluid_par.in
    filen='./input_template/fluid_par.in';
    fileID=fopen(filen);        
    FC=textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);

    str=sprintf('%05.3f',l3e);
    FC{1}{2}(1:length(str))=str;
    str=sprintf('%05.3f',l2e);
    FC{1}{3}(1:length(str))=str;
    str=sprintf('%i, %i',n3,n2);
    FC{1}{4}(1:length(str))=str;
    fID=fopen(strcat(infold,'fluid_par.in'),'w');
    for ii=1:length(FC{1})
        fprintf(fID,'%s \n',FC{1}{ii});
    end
    fclose(fID);


    % create input files form template: global_par.in
    filen='./input_template/global_par.in';
    fileID=fopen(filen);        
    FC=textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);

    fID=fopen(strcat(infold,'global_par.in'),'w');
    for ii=1:length(FC{1})
        fprintf(fID,'%s \n',FC{1}{ii});
    end
    fclose(fID);


    % copy submission script
    dest=strcat('../',cname,'/submit.sh');
    status=copyfile('./input_template/submit.sh',dest);    
    
end


