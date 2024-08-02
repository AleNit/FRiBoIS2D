
% read center of mass kinematics

function [t,p,v,a]=getcm(filename)

fileID = fopen(filename,'r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f';
size = [10 Inf];
C = fscanf(fileID,formatSpec,size)';
fclose(fileID);

t=C(:,1);
p=squeeze(C(:,2:4));
v=squeeze(C(:,5:7));
a=squeeze(C(:,8:10));


end
