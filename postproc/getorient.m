
% create orientation arrows

function [o1,o2,v1,v2]=getorient(cmpos)

o1(1)=cmpos(1);
o1(2)=cmpos(2);
v1(1)=cmpos(1);
v1(2)=cmpos(2);

r=0.5;

%apply rotation
a=cmpos(3);
v2=[v1(1)+r*sin(a),v1(2)+r*cos(a)];
o2=[o1(1)+r*cos(a),o1(2)-r*sin(a)];

end