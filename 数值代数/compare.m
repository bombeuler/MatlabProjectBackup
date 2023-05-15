function err=compare(y1,y)
n=length(y);
m=zeros(n,1);
for i=1:n
    m(i)=abs(y1(i)-y(i));
end
err=max(m);
end