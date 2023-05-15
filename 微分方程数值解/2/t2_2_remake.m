% Matlab R2022a

N=40;
h=1/N;

I1=-eye(N-1);
A1=gallery('tridiag',N-1,-1,2,-1);
Am=gallery('tridiag',N+1,-1,2,-1);
M1=kron(eye(N+1),A1)+kron(Am,eye(N-1));
M2=gallery('poisson',2*N-1);

%% M
M=zeros((3*N-1)*(N-1));

M(1:(N-1)*N,1:(N-1)*(N+1))=M1(1:(N-1)*N,1:(N-1)*(N+1));

mid_l=N*(N-1);

M(mid_l+1:mid_l+N-1,(N-1)^2+1:N*(N-1))=I1;

M(mid_l+1:mid_l+(2*N-1)*(N-1),mid_l+1:mid_l+(2*N-1)*(N-1))=M2(1:(2*N-1)*(N-1),1:(2*N-1)*(N-1));

M=M./(h^2);

%% F

F1=ones(N-1);
Fm=ones(1,N-1);
F2=ones(N-1,2*N-1);

Fl=[reshape(F1',(N-1)^2,1) ;Fm';reshape(F2',(N-1)*(2*N-1),1)];

ul_in=M\Fl;

%% uu

u1l_in=ul_in(1:(N-1)^2);
u2l_in=ul_in((N-1)^2+1:N*(N-1))';
u3l_in=ul_in(mid_l+1:end);

uu1_in=reshape(u1l_in,N-1,N-1)';
uu2_in=u2l_in;
uu3_in=reshape(u3l_in,2*N-1,N-1)';

uu1=zeros(N+1,N+1);
uu1(2:N,2:N)=uu1_in;
uu1(N+1,:)=[0,uu2_in,0];

uu2=zeros(1,2*N+1);
uu2(2:N)=uu2_in;

uu3=zeros(N,2*N+1);
uu3(1:N-1,2:2*N)=uu3_in;

uu23=[uu2;uu3];

%% mesh
xx1=-1:h:0;
yy1=-1:h:0;
xx23=-1:h:1;
yy23=0:h:1;

[X1,Y1]=meshgrid(xx1,yy1);
[X23,Y23]=meshgrid(xx23,yy23);

figure(1)
mesh(X1,Y1,uu1);
hold on
mesh(X23,Y23,uu23);
hold off
% mesh([X1,0;X23],[Y1,0;Y23],[uu1,0;uu23]);



