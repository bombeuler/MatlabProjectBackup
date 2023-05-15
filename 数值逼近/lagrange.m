function p=lagrange(x,y,t)
n=length(x);
p=0;
for i=1:n
      l=1;
    for j=1:n
        if j~=i
            l=l.*(t-x(j))./(x(i)-x(j));
        end
    end
    p=p+y(i).*l;
end

end