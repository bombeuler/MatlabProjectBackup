function lf = P2IT12(ff)
xx = -1:0.01:1;
M = diag([2,2/3,2/5]);
b = zeros(3,1);
b(1) = integral(@(x) ff(x),-1,1);
b(2) = integral(@(x) ff(x).*x,-1,1);
b(3) = integral(@(x) ff(x).*(3.*x -1)./2,-1,1);
a = M\b;
disp(a);
% [1,x,(3.*x.^2 - 1)./2]
lf = @(x) a(1) + a(2).*x +a(3).*(3.*x.^2 - 1)./2;
plot(xx,ff(xx),xx,lf(xx));
hold on
end