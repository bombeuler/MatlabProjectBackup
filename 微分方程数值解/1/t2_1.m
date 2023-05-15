%t2_1.m

%funs
f=@(x,y) -2.*exp(x+y);
g=@(x,y) exp(x+y);

%params
N=10;
Nx=N;
Ny=N;

%primary
hx=1/Nx;
hy=1/Ny;
xx=0:hx:1;
yy=0:hy:1;
[X,Y]=meshgrid(xx,yy);
uu=zeros(Nx+1,Ny+1);
uu_r=exp(X+Y);

%uu border
for i=1:Nx+1
    uu(i,1)=g(xx(i),0);
    uu(i,Ny+1)=g(xx(i),1);
end

for j=2:Ny
    uu(Nx+1,j)=g(1,yy(j));
    uu(1,j)=g(0,yy(j));
end

A=gallery('tridiag',Nx-1,-1,2,-1);
B=gallery('tridiag',Ny-1,-1,2,-1);
Ix=eye(Nx-1);Iy=eye(Ny-1);
M=kron(B./(hy^2),Ix)+kron(Iy,A./(hx^2));

%处理F数组
%
% F1,1          F1,~            F1,Ny-1
% 
%F~,1           F_in            F~,Ny-1
% 
%FNx-1,1     FNx-1,~     FNx-1,Ny-1

F=f(X(2:Nx,2:Ny),Y(2:Nx,2:Ny));

F(1,1)=F(1,1)+uu(1,2)/(hx^2)+uu(2,1)/(hy^2);
F(1,Ny-1)=F(1,Ny-1)+uu(1,Ny)/(hx^2)+uu(2,Ny+1)/(hy^2);
F(Nx-1,1)=F(Nx-1,1)+uu(Nx+1,2)/(hx^2)+uu(Nx,1)/(hy^2);
F(Nx-1,Ny-1)=F(Nx-1,Ny-1)+uu(Nx+1,Ny)/(hx^2)+uu(Nx,Ny+1)/(hy^2);

for j=2:Ny-2
    F(1,j)=F(1,j)+uu(1,j+1)/(hx^2);
    F(Nx-1,j)=F(Nx-1,j)+uu(Nx+1,j+1)/(hx^2);
end

for i=2:Nx-2
    F(i,1)=F(i,1)+uu(i+1,1)/(hy^2);
    F(i,Ny-1)=F(i,Ny-1)+uu(i+1,Ny+1)/(hy^2);
end

Fl=reshape(F,(Nx-1)*(Ny-1),1);

ul=M\Fl;
uu(2:Nx,2:Ny)=reshape(ul,Ny-1,Nx-1);

disp(disperse_norm2(uu_r-uu,2));

