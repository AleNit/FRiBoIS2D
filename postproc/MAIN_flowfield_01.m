
% read simulation output and plot FSI frames

% Alessandro Nitti, Polytechnic University of Bari (2024)

clc
clear 
close all
clearAllMemoizedCaches
kf=1;


%% input parameters
fold='../test_GallRect1DOF/out_01/';        % test folder
% fold='../test_OscCyl1DOF/out_01/';
ni=1;                                       % start time step index
ne=2000;                                    % end time step index
var=4;      % flow field: 1:z-velocity, 2:y-velocity, 3:pressure, 4:vorticity
wrf=false(1);                                % write plot to file


%% read fluid grid data
fgn=strcat(fold,'movie/nodes_grid.h5');
fgc=strcat(fold,'movie/cell_center_grid.h5');
[n2,n3,yn,zn]=readgrid_n(fgn);
[n2m,n3m,yc,zc]=readgrid_c(fgc);
zlim=zn(end);
ylim=yn(end);

[Zn,Yn]=meshgrid(zn,yn);
[Zc,Yc]=meshgrid(zc,yc);


%% read initial body geometry
[zp,yp,zcm0,ycm0]=getigeo(fold);


%% read body kinematics
[t,pcm,vcm,acm]=getcm(strcat(fold,'movie/movierb.out'));


%% start time loop
figure(kf); kf=kf+1; 
% figure('visible','off');
set(gcf,'Position',[100 100 700 700])
T=tiledlayout(3,1);

if (wrf)
    if (~exist('./frames/', 'dir'))
       mkdir('./frames/')
    end
end

for n=ni:1:ne
    
    fframe=strcat(fold,'movie/frame_',sprintf('%05d',n),'.h5');
    
    [vn,wn,pc]=getflowdata(fframe);
    
    [vrtx]=vorticity2d(zc,yc,wn',vn',n3,n2);
    
    % plot contours  
    nexttile(1,[2,1])
    switch var
        case 1
            vart='z-velocity';
            surf(Zn(1:end-1,1:end-1),Yc,zeros(size(Yc)),wn(1:end-1,1:end-1));   
            colormap(parula); clim([-2,2])
        case 2
            vart='y-velocity';
            surf(Zc,Yn(1:end-1,1:end-1),zeros(size(Zc)),vn(1:end-1,1:end-1));        
            colormap(parula); clim([-2,2])
        case 3
            vart='pressure';
            surf(Zc,Yc,zeros(size(Zc)),pc);        
            colormap(parula); clim([-5,5])
        case 4
            vart='vorticity';
            surf(Zc,Yc,zeros(size(Zc)),vrtx');             
            colormap(getmap); clim([-3,3])
    end
    hold on    
    shading interp  
    colorbar    
    
    % plot filled geometry
    delete(findobj('type', 'patch'));
    delete(findobj('type', 'annotation'));
    zcm=pcm(n,2);
    ycm=pcm(n,1);
    ang=pcm(n,3);
    zpr=(zp-zcm0).*cos(-ang) - (yp-ycm0).*sin(-ang) + zcm;
    ypr=(zp-zcm0).*sin(-ang) + (yp-ycm0).*cos(-ang) + ycm;
    patch(zpr,ypr,'white')
    
    % % plot orientation axes
    % [o1,o2,v1,v2]=getorient([zcm,ycm,ang]);
    % p1=[o1(1) o1(2) 0];
    % p2=[o2(1) o2(2) 0];
    % arrow3d(p1,p2,20,'line',[0.5,0.5],[10,10],'k');
    % p1=[v1(1) v1(2) 0];
    % p2=[v2(1) v2(2) 0];
    % arrow3d(p1,p2,20,'line',[0.5,0.5],[10,10],'k');
    
    axis on    
    axis equal
    % axis([zcm0-4,zcm0+10,ycm0-4,ycm0+4])          % close-up (cylinder)
    axis([zcm0-7,zcm0+14,ycm0-7,ycm0+7])          % close-up (rectangle)
    % axis([0,zlim,0,ylim])        % full domain
    xlabel('$z/D$','interpreter','latex','fontsize',14);
    ylabel('$y/D$','interpreter','latex','fontsize',14);  
    strt=strcat(vart,', tU/D=',sprintf('%4.2f',t(n)));
    title(strt,'interpreter','latex','fontsize',14)
    box on
    set(gca,'BoxStyle','full')
    view(2)
    lighting none

    drawnow
    
    % draw verification line values
    nexttile(3,[1,1])        
    if (isequal(fold,'../test_OscCyl1DOF/out_01/'))
        plot(t(1:n),pcm(1:n,1)-30,'-b')
        hold on
        plot(t(1:n),0.57.*ones(n,1),'--r')
        plot(t(1:n),-0.57.*ones(n,1),'--r','HandleVisibility','off')        
        ylabel('$(y-y_0)/D$','interpreter','latex','fontsize',14); 
        legend('present','Bourguet, LoJacono (2014)','interpreter','latex', ...
            'fontsize',12,'orientation','horizontal'); 
        axis([0,t(n)+2,-1.5,1.5])
    elseif (isequal(fold,'../test_GallRect1DOF/out_01/'))
        plot(t(1:n),pcm(1:n,3).*180/pi,'-b')
        hold on
        plot(t(1:n),12.2.*ones(n,1),'--r')
        plot(t(1:n),-12.2.*ones(n,1),'--r','HandleVisibility','off')                
        ylabel('$\theta$','interpreter','latex','fontsize',14);
        legend('present','Robertson (2002)','interpreter','latex', ...
            'fontsize',12,'orientation','horizontal'); 
        axis([0,t(n)+2,-25,25])
    end
    xlabel('$t \, U/D$','interpreter','latex','fontsize',14);
        
    drawnow
    hold off


    if (wrf)
        picname=strcat('./frames/snap_',sprintf('%05d',n));
        print('-dpng','-r200',picname);                        
    end
    
    hold off
    disp(['step ' num2str(n)])

 end


