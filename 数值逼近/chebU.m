%�ڶ����б�ѩ�����ʽ�ƽ�
function p=chebU(f,n,s)
p=0;
for j=0:n
    p=p+(2/pi).*integral(@(x) sin((j+1).*acos(x)).*f(x),-1,1).*sin((j+1).*acos(s))./sqrt(1-s.^2);
end
end