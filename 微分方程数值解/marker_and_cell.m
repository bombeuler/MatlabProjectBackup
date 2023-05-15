h=0.05;
N=1/h;

A1=gallery('tridiag',N,-1,2,-1);
A2=gallery('tridiag',N-1,-1,2,-1);
A3=gallery('tridiag',N,-1,2,-1);
A3(1,1)=3;
A3(N,N)=3;
A4=gallery('tridiag',N-1,-1,2,-1);
A4(1,1)=3;
A4(N-1,N-1)=3;
% A4=blkdiag(1,A4,1);
I1=speye(N-1);
I2=speye(N);
I3=speye(N+1);
O1=zeros(N*(N-1));
O2=zeros(N*N,N*N);
% I3(1,1)=0;
% I3(N+1,N+1)=0;

M1=(kron(I1,A3)+kron(A2,I2));
M2=(kron(I2,A4)+kron(A1,I1));

% M1=blkdiag(I2,M1,I2);

P1=gallery('tridiag',N,0,-1,1);
Q1=kron(P1,I2);
Q1=Q1(1:N*(N-1),:);

P2=gallery('tridiag',N,0,-1,1);
P2=P2(1:end-1,:);
Q2=kron(I2,P2);

P3=gallery('tridiag',N,-1,1,0);
U1=kron(P3,I2);
U1=U1(:,1:(N-1)*N);

P4=-P2';
U2=kron(I2,P4);

M1=M1./(h^2);
M2=M2./(h^2);
Q1=Q1./h;
Q2=Q2./h;
U1=U1./h;
U2=U2./h;

F1=zeros(N,N-1);
F1(N,:)=2./(h^2);
F1l=reshape(F1,N*(N-1),1);
F2=zeros(N-1,N);
F2l=reshape(F2,N*(N-1),1);
F3l=zeros(N*N,1);

M=[M1,O1,Q1;O1,M2,Q2;U1,U2,O2];
M(2*N*(N-1)+1,:)=0;
M(2*N*(N-1)+1,2*N*(N-1)+1)=1;
Fl=[F1l;F2l;F3l];

upl=M\Fl;

%Jacobi
% D=diag(diag(M));
% 
% tol=0.01;
% 
% u_old=ones(length(diag(M)),1);
% u_new=-inv(D)*(M-D)*u_old+D\Fl;
% disp(u_new);
% 
% while (norm(u_new-u_old,2)>=tol)
%     u_old=u_new;
%     disp('jj');
%     u_new=-inv(D)*(M-D)*u_old+D\Fl;
% end
% 
% upl=u_new;

u1l=upl(1:N*(N-1));
u2l=upl(N*(N-1)+1:2*N*(N-1));
pl=upl(2*N*(N-1)+1:end);

uu1=zeros(N,N+1);
uu2=zeros(N+1,N);
uu1(:,2:N)=reshape(u1l,N,N-1);
uu2(2:N,:)=reshape(u2l,N-1,N);

pp=zeros(N,N);
pp(1:N,1:N)=reshape(pl,N,N);

uu1m=zeros(N);
uu2m=zeros(N);

for i=1:N
    uu1m(:,i)=(uu1(:,i)+uu1(:,i+1))./2;
    uu2m(i,:)=(uu2(i,:)+uu2(i+1,:))./2;
end

u1=zeros(N+2);
u2=zeros(N+2);

u1(N+2,:)=1;
u1(2:N+1,2:N+1)=uu1m;
u2(2:N+1,2:N+1)=uu2m;

xx1=[0,h/2:h:1-h/2,1];
yy1=[0,h/2:h:1-h/2,1];

[X,Y]=meshgrid(xx1,yy1);

figure(1)
quiver(X,Y,u1,u2);

% xx1=0:h:1;
% yy1=[0,h/2:h:1-h/2,1];
% 
% [X1,Y1]=meshgrid(xx1,yy1);
% 
% xx2=[0,h/2:h:1-h/2,1];
% yy2=0:h:1;
% 
% [X2,Y2]=meshgrid(xx2,yy2);
% 
xx3=h/2:h:1-h/2;
yy3=h/2:h:1-h/2;

[X3,Y3]=meshgrid(xx3,yy3);

% figure(1);
% mesh(X1,Y1,uu1);
% 
% figure(2);
% mesh(X2,Y2,uu2);
% 
figure(2);
mesh(X3,Y3,pp);

