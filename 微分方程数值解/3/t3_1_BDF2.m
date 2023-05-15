N=10;
Tm=1;
h=1/N;
dt=0.0001;

f=@(x,y,t) (1+2.*pi.^2).*exp(t).*sin(pi.*x).*sin(pi.*y);
g=@(x,y,t) x.*0+y.*0+t.*0;
u0=@(x,y) sin(pi.*x).*sin(pi.*y);
ur_fun=@(x,y,t) exp(t).*sin(pi.*x).*sin(pi.*y);

A=gallery('poisson',N-1);
I=speye((N-1)*(N-1));

M1=2*I+A.*dt./(h^2);
M2=3/2*I-A.*dt./(h^2);

xx=0:h:1;
yy=0:h:1;
tt=0:dt:Tm;
[X,Y,T]=meshgrid(xx,yy,tt);
U=g(X,Y,T);
Ur=ur_fun(X,Y,T);

F=f(X,Y,T);

u_older=zeros((N-1)*(N-1),1);
u_old=zeros((N-1)*(N-1),1);
u_new=zeros((N-1)*(N-1),1);

for ti=1:length(tt)
    if ti==1
       U(:,:,ti)=u0(X(:,:,1),Y(:,:,1));
       u_older=reshape(U(2:N,2:N,ti)',(N-1)*(N-1),1);
    elseif ti==2
        uu=U(:,:,ti);
        F_square=F(2:N,2:N,ti);
        F_square(1,1)=F_square(1,1)+(uu(1,2)+uu(2,1))./(h^2);
        F_square(1,N-1)=F_square(1,N-1)+(uu(2,N+1)+uu(1,N))./(h^2);
        F_square(N-1,1)=F_square(N-1,1)+(uu(N+1,2)+uu(N,1))./(h^2);
        F_square(N-1,N-1)=F_square(N-1,N-1)+(uu(N,N+1)+uu(N+1,N))./(h^2);
        F_square(1,2:N-2)=F_square(1,2:N-2)+uu(1,3:N-1)./(h^2);
        F_square(N-1,2:N-2)=F_square(N-1,2:N-2)+uu(N+1,3:N-1)./(h^2);
        F_square(2:N-2,1)=F_square(2:N-2,1)+uu(3:N-1,1)/(h^2);
        F_square(2:N-2,N-1)=F_square(2:N-2,N-1)+uu(3:N-1,N+1)./(h^2);

        Fl=reshape(F_square',(N-1)*(N-1),1);
        u_old=M1*u_older+dt.*Fl;
        U(2:N,2:N,ti)=reshape(u_new,N-1,N-1)';
    else
         uu=U(:,:,ti);
        F_square=F(2:N,2:N,ti);
        F_square(1,1)=F_square(1,1)+(uu(1,2)+uu(2,1))./(h^2);
        F_square(1,N-1)=F_square(1,N-1)+(uu(2,N+1)+uu(1,N))./(h^2);
        F_square(N-1,1)=F_square(N-1,1)+(uu(N+1,2)+uu(N,1))./(h^2);
        F_square(N-1,N-1)=F_square(N-1,N-1)+(uu(N,N+1)+uu(N+1,N))./(h^2);
        F_square(1,2:N-2)=F_square(1,2:N-2)+uu(1,3:N-1)./(h^2);
        F_square(N-1,2:N-2)=F_square(N-1,2:N-2)+uu(N+1,3:N-1)./(h^2);
        F_square(2:N-2,1)=F_square(2:N-2,1)+uu(3:N-1,1)/(h^2);
        F_square(2:N-2,N-1)=F_square(2:N-2,N-1)+uu(3:N-1,N+1)./(h^2);

        Fl=reshape(F_square',(N-1)*(N-1),1);
        
        u_new=M2\(2.*u_old-u_older./2);
        
        u_older=u_old;
        u_old=u_new;
        U(2:N,2:N,ti)=reshape(u_new,N-1,N-1)';
    end
end

disp(max(max(max(abs(U-Ur)))));
