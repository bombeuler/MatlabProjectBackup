function [x,Q,R]=zqrfact(A,b)
[m,n]=size(A);
Q=eye(m);
R=zeros(n);
for j=1:n
    [v,beta]=householder(A(j:m,j));
    A(j:m,j:n)=(eye(m-j+1)-beta*v*v')*A(j:m,j:n);
    d(j)=beta;
    if j<m
        A(j+1:m,j)=v(2:m-j+1);
    end
end
R=triu(A);
for k=1:n
    H=eye(m);
    H1=eye(m-k+1)-d(k)*[1,A(k+1:m,k)']'*[1,A(k+1:m,k)'];
    H(k:m,k:m)=H1;
    Q=Q*H;
end
x=zeros(n,1);
b=Q'*b;
for j=n:-1:1
    x(j)=(b(j)-sum(R(j,j+1:n).*x(j+1:n)'))/R(j,j);
end
end