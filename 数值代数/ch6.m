% [k,c]=half(3,4,10^(-8));
% disp(c);
format long;

% disp(k);
% function [k,c]=half(a,b,tol)
% c=(a+b)/2;k=1;
% m=1+round((log(b-a)-log(2*tol))/log(2));
% while k<=m
%     if f(c)==0
%         c=c;
%         return;
%     elseif f(a)*f(c)<0
%             b=(a+b)/2;
%         else
%             a=(a+b)/2;
%         end
%         c=(a+b)/2;
%         k=k+1;
%     end
% 
% end
% function y=f(z)
%  y=z^3-48;
% end
% [k,c]=steffensen(pi/4,0,10^(-8));
% disp(k);
% disp(c);
% d=compare([1;2;3;4],[0;4;2;1]);
% disp(d)
x=[7*pi/6,4*pi/3;];
disp(x-2*sin(x+pi/3));
disp(newton(pi/2,10^(-8)));
function [k,c]=secant(a,b,tol)
k=0;
while abs(b-a)>tol
c=b-(b-a)*f(b)/(f(b)-f(a));
a=b;
b=c;
k=k+1;
end
end1
function [k,c]=steffensen(a,b,tol)
k=0;
while abs(b-a)>tol
    c=a-f(a)^2/(f(a)-f(a-f(a)));
    b=a;
    a=c;
    k=k+1;
end
end
function y=f(x)
y=2*cos(x)-1-sin(x);
end
function x=picard(x0,tol)
    x1=y(x0);
    while(abs(x1-x0)>=tol)
        x0=x1;
        x1=y(x0);
    end
    x=x1;
end
function m=y(x)
a=cos(1/2);
    m=(a*x-sin(x)+1)/(1+a);
end
function x=newton(x0,tol)
    x1=nf(x0);
    k=1;
    while(abs(x1-x0)>=tol)
        k=k+1;
        x0=x1;
        x1=nf(x0);
    end
    x=x1;
end
function y=nf(x)
    y=x-(x-2*sin(x+pi/3))/(1-2*cos(x+pi/3));
end