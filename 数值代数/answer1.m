format short;
%%确定A,b,h
h=0.1;
m1=ones(8,1);
m2=ones(9,1).*(4*h^2-2);
m3=ones(8,1);
A=diag(m2)+diag(m1,-1)+diag(m3,1);
b=zeros(9,1);
b(1)=r(0.1)*h^2-2;
b(9)=r(0.9)*h^2+2;
for i=2:8
    b(i)=r(i*0.1)*h^2;
end
disp(A);
disp(b);
%原值
fl=zeros(9,1);
for i=1:9
    fl(i)=f(0.1*i);
end
disp(fl);

%Doolitle
[y1,L,U]=lux(A,b);
disp(y1);
err1 =compare(y1,fl);
disp(err1); 

%Cholesky
y2=cholesky(A,b);
disp(y2);
err2=compare(y2,fl);
disp(err2); 

%QR
[y3,Q,R]=qrfact(A,b);
disp(y3);
err3 =compare(y3,fl);
disp(err3); 

%追赶法
y4=chase(m1,m2,m3,b);
disp(y4);
err4 =compare(y4,fl);
disp(err4); 

function c=r(x)
c=(4-pi^2)*(2*cos(pi*x)+3*sin(pi*x));
end
function y=f(x)
y=2*cos(pi*x)+3*sin(pi*x);
end
function [v,beta]= householder(x)
n=length(x);eta=norm(x,inf);x=x/eta;
sigma=x(2:n)'*x(2:n); v=x; v(1)=1;
if sigma==0
    beta=0;
else 
    alpha=sqrt(x(1)^2+sigma);
    if x(1)<=0
        v(1)=x(1)-alpha;
    else
        v(1)=-sigma/(x(1)+alpha);
    end
    beta=2*v(1)^2/(sigma+v(1)^2);
    v=v/v(1);
end
end
%% Doolitle
function [x,L,U]=lux(A,b)
[n,n]=size(A);
L=zeros(n);
U=zeros(n);
x=zeros(n,1);
y=zeros(n,1);
for r=1:n
    for i=r:n
        U(r,i)=A(r,i)-sum(L(r,1:r-1).*U(1:r-1,i)');
        L(i,r)=(A(i,r)-sum(L(i,1:r-1).*U(1:r-1,r)'))/U(r,r);
    end
end
for i=1:n
    y(i)=b(i)-sum(L(i,1:i-1).*y(1:i-1)');
end
for j=n:-1:1
    x(j)=(y(j)-sum(U(j,j+1:n).*x(j+1:n)'))/U(j,j);
end
end
%% Cholesky
function x=cholesky(a,b)
n=length(b);
v=zeros(n);
x=zeros(n,1);
y=zeros(n,1);
for j=1:n
    for i=1:j-1
        v(j,i)=a(j,i)*a(i,i);
    end
    a(j,j)=a(j,j)-a(j,1:j-1)*v(j,1:j-1)';
    a(j+1:n,j)=(a(j+1:n,j)-a(j+1:n,1:j-1)*v(j,1:j-1)')/a(j,j);
end
L=tril(a,-1)+eye(n);
U=diag(diag(a))*L';
for i=1:n
    y(i)=b(i)-L(i,1:i-1)*y(1:i-1);
end
for j=n:-1:1
    x(j)=(y(j)-U(j,j+1:n)*x(j+1:n))/U(j,j);
end
end
%% QR
function [x,Q,R]=qrfact(A,b)
n=length(A); 
Q=eye(n); 
R=zeros(n);
for j=1:n
    if j<n
        [v,beta]=householder(A(j:n,j));
    else
        v=1; beta=2-2*mod(n,2);
    end
    A(j:n,j:n)=(eye(n-j+1)-beta*v*v')*A(j:n,j:n);
    d(j)=beta;
    if j<n
        A(j+1:n,j)=v(2:n-j+1);
    end
end
R=triu(A);
for k=1:n
    H=eye(n);
    H1=eye(n-k+1)-d(k)*[1,A(k+1:n,k)']'*[1,A(k+1:n,k)'];
    H(k:n,k:n)=H1;
    Q=Q*H;
end
x=zeros(n,1);b=Q'*b;
for j=n:-1:1
    x(j)=(b(j)-sum(R(j,j+1:n).*x(j+1:n)'))/R(j,j);
end
end
%% 追赶法
function x=chase(a,b,c,d)
n=length(b);
f(1)=c(1)/b(1);
g(1)=d(1)/b(1);
for i=2:n-1
    h(i)=b(i)-f(i-1)*a(i-1);
    f(i)=c(i)/h(i);
    g(i)=(d(i)-g(i-1)*a(i-1))/h(i);
end
g(n)=(d(n)-g(n-1)*a(n-1))/(b(n)-f(n-1)*a(n-1));
x(n)=g(n);
for i=n-1:-1:1
    x(i)=g(i)-f(i)*x(i+1);
end
x=x';
end
%% 误差函数
function err=compare(y1,y)
n=length(y);
m=zeros(n,1);
for i=1:n
    m(i)=abs(y1(i)-y(i));
end
err=max(m);
end