
function xn=cubicstr(int,nint,li,le)

int=int-li;
le=le-li;

% uniform spacing region
xnm=linspace(int(1),int(2),nint(2)+2);
xnm(1)=[];
xnm(end)=[];

% stretching region 2
dx0=(int(2)-int(1))/(nint(2)+2);
B=dx0*nint(3);
A=(le-int(2))-B;
xn2=zeros(1,nint(3));
s=linspace(0,1,nint(3));
for i=1:nint(3)
    xn2(i)=A*s(i)^3+B*s(i);
end
xn2=xn2+int(2);

% stretching region 1
B=dx0*nint(1);
A=int(1)-B;
xn1=zeros(1,nint(1));
s=linspace(0,1,nint(1));
for i=1:nint(1)    
    xn1(i)=A*s(i)^3+B*s(i);
end
xn1=-xn1+int(1);
xn1=flip(xn1);

xn=[xn1 xnm xn2];
xn=xn+li;


end
