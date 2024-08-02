
function [v,w,pr]=getflowdata(filename)

fileID=strcat(filename);

fid = H5F.open(fileID);

dset_id = H5D.open(fid,'Vy');
v = H5D.read(dset_id);
H5D.close(dset_id);

dset_id = H5D.open(fid,'Vz');
w = H5D.read(dset_id);
H5D.close(dset_id);

dset_id = H5D.open(fid,'Pr');
pr = H5D.read(dset_id);
H5D.close(dset_id);

H5F.close(fid);

end

