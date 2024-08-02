
function myMap=getmap

c_bot = [0,0.4470,0.7410];
c_mid = [1.0,1.0,1.0];
c_top = [0.850,0.3250,0.098];
mp1 = [linspace(c_bot(1),c_mid(1),50)',linspace(c_bot(2),c_mid(2),50)',linspace(c_bot(3),c_mid(3),50)'];
mp2 = [linspace(c_mid(1),c_top(1),50)',linspace(c_mid(2),c_top(2),50)',linspace(c_mid(3),c_top(3),50)'];
myMap=[mp1;mp2];

end