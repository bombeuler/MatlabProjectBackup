function plot2LambdaT1(x0,x1)
lambda0 = @(x) (x1-x)./(x1-x0);
lambda1 = @(x) (x-x0)./(x1-x0);
xx = x0:0.01:x1;
figure(1);
hold on
plot(xx,lambda0(xx),'-b');
plot(xx,lambda1(xx),'--b');
plot(xx,lambda0(xx)+lambda1(xx),'.b');
plot(xx,x0.*lambda0(xx),'-r');
plot(xx,x1.*lambda1(xx),'--r');
plot(xx,x0.*lambda0(xx)+x1.*lambda1(xx),'.r');
end