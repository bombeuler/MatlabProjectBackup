function NewtonBiStableSolver2D()
g = [2     2     2     2
     0     1     1     0
     1     1     0     0
     0     0     1     1
     0     1     1     0
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'hmax',0.02);
x=p(1,:); y=p(2,:);
xi_old=(cos(2*pi*x.^2).*cos(2*pi*y.^2))';
xi_new=xi_old;
dt=0.1; % time step
epsilon=0.01;
[A,M]=assema(p,t,1,1,0);

v = VideoWriter('bistable');
% v.FileFormat='mp4';
v.FrameRate=100;
open(v);

for l=1:300 % time loop
for k=1:3 % non-linear loop
ii=t(1,:); jj=t(2,:); kk=t(3,:);
xi_tmp=xi_new; % copy temporary solution to new
xi_tmp_mid=(xi_tmp(ii)+xi_tmp(jj)+xi_tmp(kk))/3;
f =(xi_tmp_mid-xi_tmp_mid.^3); % evaluate f
df=1-3*xi_tmp_mid.^2; % evaluate derivative df of f
[crap,Mdf,b]=assema(p,t,0,df',f');
J=(M+dt*epsilon*A)-dt*Mdf; % Jacobian
rho=(M+dt*epsilon*A)*xi_new ...
-M*xi_old-dt*b; % residual
xi_new=xi_tmp-J\rho; % Newton update
error=norm(xi_tmp-xi_new);
end
xi_old=xi_new; % copy old solution to new
pdesurf(p,t,xi_new)
axis([0 1 0 1 -1 1]), caxis([-1,1]);
frame = getframe(gcf);
writeVideo(v,frame);
end


close(v);