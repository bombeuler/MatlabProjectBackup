ee=1e-5;

u=@(x) x-(exp(-(1-x)/ee)-exp(-1/ee))/(1-exp(-1/ee));

for i=0:22
N=10*(2^i);
h=1/N;
uu=zeros(1,N+1);
% A=gallery('tridiag',N-1,-ee/(h^2)-1/(2*h),2*ee/(h^2),1/(2*h)-ee/(h^2));
A=gallery('tridiag',N-1,-1/h-ee/(h^2),1/h+2*ee/(h^2),-ee/(h^2));
Fl=ones(N-1,1);
uu(2:N)=(A\Fl)';

xx=0:h:1;
ur=u(xx);
err_l=abs(ur-uu);
err=max(err_l);
% disp(err);
plot(i,err,'or');
hold on
% figure(1);
% plot(xx,ur);
% figure(2);
% plot(xx,err_l,'-b');
end
hold off