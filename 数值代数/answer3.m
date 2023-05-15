format short;

%% 初始化A,b
A=zeros(50,50);
b=zeros(50,1);
for i=1:50
    for j=1:50
        if(i==j)
            A(i,j)=50*i;
        else
            A(i,j)=max([i,j]);
        end
    end
end
for l=1:50
    for k=1:50
        b(l)=b(l)+A(l,k)*(51-k);
    end
end
disp(A);
disp(b);
M=diag(diag(A));
%% 进行计算
X_1=sdescent(A,b);
disp(X_1);
X_2=congrad(A,b);
disp(X_2);
X_3=precongrad(A,b,M);
disp(X_3);

%%  不同方法函数

%最速下降法
function X=sdescent(A,b)
n=length(b);
X0=zeros(n,1);
r0=b-A*X0;
k=0;
while norm(r0,inf)>=10^(-8)
    k=k+1;
    alpha0=r0'*r0/(r0'*A*r0);
    X1=X0+alpha0*r0;
    r0=b-A*X1;
    X0=X1;
end
X=X0;
end

%共轭梯度法
function X=congrad(A,b)
n=length(b);
X0=zeros(n,1);
r0=b-A*X0;
P0=r0;
while norm(r0,inf)>=10^(-8)
    alpha=r0'*r0/(P0'*A*P0);
    X1=X0+alpha*P0;
    r1=r0-alpha*A*P0;
    beta=r1'*r1/(r0'*r0);
    P1=r1+beta*P0;
    r0=r1;
    X0=X1;
    P0=P1;
end
X=X0;
end

%预优共轭梯度法
function X=precongrad(A,b,M)
n=length(b);
X0=zeros(n,1);
r0=b-A*X0;
zeta0=M\r0;
rho0=r0'*zeta0;
P0=zeta0;
while norm(r0,inf)>=10^(-8)
    omega=A*P0;
    alpha=rho0/(P0'*omega);
    X1=X0+alpha*P0;
    r1=r0-alpha*omega;
    zeta1=M\r1;
    rho1=r1'*zeta1;
    beta=rho1/rho0;
    P1=zeta1+beta*P0;
    r0=r1;
   X0=X1;
   P0=P1;
   rho0=rho1;
end
X=X0;
end


