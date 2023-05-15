g = Dslitg();
h=0.025;
k=0.005;
T=5;
[p,e,t] = initmesh(g,"hmax",h);
np = size(p,2);
x = p(1,:)';
y = p(2,:)';
fixed = find(x< -0.24999);
xi = zeros(np,1);
eta = zeros(np,1);

[A,M,b] = assema(p,t,1,1,0);

v = VideoWriter('peaks.mp4');
open(v);

for l = 1:round(T/k)
    time = l*k;
    LHS = [M -0.5*k*M; 0.5*k*A M];
    rhs = [M 0.5*k*M; -0.5*k*A M]*[xi; eta] +[zeros(np,1); k*b];
    sol = LHS\rhs;
    xi = sol(1:np);
    eta = sol(np+1:end);
    xi(fixed) = 0.1 * sin(8*pi*time);
    pdeplot(p,[],t,"XYData",xi,"ZData",xi);
    axis([-0.3 1 0 1 -0.25 0.5]);
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
