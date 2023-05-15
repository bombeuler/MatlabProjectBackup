function p = newton(x,y,t)
n=length(x);
F=zeros(n,n);
F(1,:)=y;
for i=2:n
    for j=2:n
        if j>=i
            F(i,j)=(F(i-1,j-1)-F(i-1,j))/(x(j-i+1)-x(j));
        end
    end
end
a=diag(F);
p=a(n);
for m=n-1:-1:1
    p=a(m)+p.*(t-x(m));
end
end