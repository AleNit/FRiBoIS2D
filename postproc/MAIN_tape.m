
% get tape from frames

clc
clear 
close all


%% video setting
outfold='./frames/';
% tapename='./OscCyl1DOF';
tapename='./GallRect1DOF';
outputVideo = VideoWriter(tapename,'MPEG-4');
outputVideo.FrameRate=35;      
outputVideo.Quality=50;         % Default 75


%% start assembly
disp('Time for animation assembly')

tic
open(outputVideo)

cd (outfold)

imagenames = dir('*.png');
imagenames = {imagenames.name}';
ax=gca();
ll=length(imagenames);

for ii=1:ll
    img=imread(fullfile(imagenames{ii}));
    writeVideo(outputVideo,img)
end

close(outputVideo)
toc

close all
