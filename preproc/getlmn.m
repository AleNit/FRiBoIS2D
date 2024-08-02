
function [zlm,ylm,Alm,norlm2]=getlmn(nlm,zcp,ycp,or,dzmin,dymin)

zlm=zeros(nlm,1);
ylm=zeros(nlm,1);
Alm=zeros(nlm,1);
norlm=zeros(nlm,2);
for i=1:nlm
    zlm(i)=(zcp(i)+zcp(i+1))*0.5;
    ylm(i)=(ycp(i)+ycp(i+1))*0.5;
    Alm(i)=sqrt((zcp(i+1)-zcp(i))^2+(ycp(i+1)-ycp(i))^2);    
end

% create normal arrays; the sign of "or" provides the orientation of the normal arrays 
dmin=min(dzmin,dymin);
for i=1:nlm
    tri_ver(1:3)=[0.0 ycp(i) zcp(i)];
    tri_ver(4:6)=[0.0 ycp(i+1) zcp(i+1)];
    tri_ver(7:9)=[or (ycp(i)+ycp(i+1))/2 (zcp(i)+zcp(i+1))/2];
    v1=tri_ver(4:6)-tri_ver(1:3);
    v2=tri_ver(7:9)-tri_ver(1:3);
    tri_nor=cross(v1,v2);
    tri_nor=tri_nor./sqrt(sum(tri_nor.^2));
    tri_nor=tri_nor.*dmin;
    norlm(i,1)=0.0+tri_nor(1);
    norlm(i,2)=ylm(i)+tri_nor(2);
    norlm(i,3)=zlm(i)+tri_nor(3);    
end

norlm2=squeeze(norlm(:,2:3));

end
