
function [n2m,n3m,yc,zc]=readgrid_c(fname)

fid = H5F.open(fname);
 
dset_id=H5D.open(fid,'rm');
yc=H5D.read(dset_id);
H5D.close(dset_id);
 
dset_id=H5D.open(fid,'zm');
zc=H5D.read(dset_id);
H5D.close(dset_id);
 
H5F.close(fid);

n2m=length(yc);
n3m=length(zc);
 

end


