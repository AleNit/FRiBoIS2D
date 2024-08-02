
function [n2,n3,yn,zn]=readgrid_n(fname)

fid = H5F.open(fname);
 
dset_id=H5D.open(fid,'rc');
yn=H5D.read(dset_id);
H5D.close(dset_id);
 
dset_id=H5D.open(fid,'zz');
zn=H5D.read(dset_id);
H5D.close(dset_id);
 
H5F.close(fid);
 
n2=length(yn);
n3=length(zn);
 

end


