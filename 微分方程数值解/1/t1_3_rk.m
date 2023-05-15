t1=10;
t2=28;
t3=8/3;
h=0.01;
T=1000;
F=@(xyz) [t1.*(xyz(2)-xyz(1)),t2.*xyz(1)-xyz(2)-xyz(1).*xyz(3),xyz(1).*xyz(2)-t3.*xyz(3)];

N=T/h;
tt=0:h:T;
uu=zeros(N+1,3);
uu(1,:)=[0.1,0,0];

for i=1:N
   ui=uu(i,:);
   k1= F(ui);
   k2=F(ui+h/2*k1);
   k3=F(ui+h/2*k1);
   uu(i+1,:)=ui+h/2*(k2+k3);
   
    
end
plot3(uu(:,1),uu(:,2),uu(:,3));