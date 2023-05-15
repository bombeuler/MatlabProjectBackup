%第一类切比雪夫多项式逼近
function p=chebT(f,n,s)
p= (integral(@(x) f(x)./sqrt(1-x.^2),-1,1)/pi);
for j=1:n
    p=p+(2/pi).*integral(@(x) f(x).*cos(j.*acos(x))./sqrt(1-x.^2),-1,1).*cos(j.*acos(s));
end
end
