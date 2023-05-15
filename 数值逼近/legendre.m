function p=legendre(f,n,s)
p=0;
for j=0:n
    p=p+((2*j+1)/2).*integral(@(x) f(x).*legendreP(j,x),-1,1).*legendreP(j,s);
end
end