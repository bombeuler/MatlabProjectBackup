%t2_1.m

%funs
f=@(x,y) -5.*exp(x-2.*y);
g=@(x,y) exp(x-2.*y);

%params
N0=10;
% Nx=10;
% Ny=20;


for fig=1:3
    N=N0^fig;
    Nx=N;
    Ny=N;
   result= t2_1_func(Nx,Ny,f,g,fig);
   figure(100);
   plot(fig,result,'or');
   hold on

end


%primary
function result=t2_1_func(Nx,Ny,f,g,fig)
hx=1/Nx;
hy=1/Ny;
xx=0:hx:1;
yy=0:hy:1;
[Y,X]=meshgrid(yy,xx);
uu=zeros(Nx+1,Ny+1);

%real
uu_r=exp(X-2.*Y);

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

%����F����
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

Fl=reshape(F',(Nx-1)*(Ny-1),1);

ul=M\Fl;
uu(2:Nx,2:Ny)=reshape(ul,Ny-1,Nx-1)';

disp(disperse_norm2(uu_r-uu,0));

figure(fig);
mesh(X,Y,uu);
result = disperse_norm2(uu_r-uu,0);
end


