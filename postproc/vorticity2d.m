
%It computes vorticity out-of-plane component from a 2D field solution
%input field are located et the cell center

function [vrtx]=vorticity2d(x,y,u,v,d1,d2)

uc=zeros(d1,d2);
vc=zeros(d1,d2);
for j=1:d2-1
   jp=j+1;   
   for i=1:d1-1
      ip=i+1;                  
      uc(i,j)=(u(ip,j)+u(i,j))*0.5;
      vc(i,j)=(v(i,jp)+v(i,j))*0.5;      
   end
end

dudy=zeros(d1-1,d2-1);
dvdx=zeros(d1-1,d2-1);
for j=1:d2-2
   jp=j+1;
   dy=1.0/(y(jp)-y(j));
   for i=1:d1-2
      ip=i+1;
      dx=1.0/(x(ip)-x(i));
            
      dudy(i,j)=(uc(i,jp)-uc(i,j))*dy;
      dvdx(i,j)=(vc(ip,j)-vc(i,j))*dx;
      
   end
end

vrtx=dudy-dvdx;

end

