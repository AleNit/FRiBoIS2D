
% read initial position of lagrangian markers and initial center of mass
% position

function [z,y,zcm0,ycm0]=getigeo(fold)

fname1=strcat(fold,'../input_FSI/lmdata.in');
fileID = fopen(fname1,'r');
formatSpec = '%f %f';
size = [2 Inf];
C = fscanf(fileID,formatSpec,size)';
fclose(fileID);
z=C(:,2);
y=C(:,1);

fname2=strcat(fold,'../input_FSI/rigid_par.in');
fileID = fopen(fname2,'r');
B=textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);
tmp=strsplit(B{1}{3});
ycm0=str2num(tmp{1});
zcm0=str2num(tmp{2});


end
