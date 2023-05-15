t1=10;
t2=28;
t3=8/3;
h=0.01;
T=1000;

N=T/h;
tt=0:h:T;
uu=zeros(N+1,3);
uu(1,:)=[0.1,0,0];

for i=1:N
    x=uu(i,1);y=uu(i,2);z=uu(i,3);
    uu(i+1,1)=x+h*t1*(y-x);
    uu(i+1,2)=y+h*(t2*x-y-x*z);
    uu(i+1,3)=z+h*(x*y-t3*z);
end
plot3(uu(:,1),uu(:,2),uu(:,3));