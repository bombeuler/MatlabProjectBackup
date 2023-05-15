% 求离散范数的函数
%0为0，1为1，c为2
function res = disperse_norm(u,flag)
N=length(u)-1;
h=1/N;
switch flag
    case 2
        res=max(abs(u));
    case 1
        res=sqrt(sum(h.*(u(2:N)).^2)+h*(u(1)^2+u(N+1)^2)/2+sum((u(2:N+1)-u(1:N)).^2./h));
    case 0
        res=sqrt(sum(h.*(u(2:N)).^2)+h*(u(1)^2+u(N+1)^2)/2);
end
end