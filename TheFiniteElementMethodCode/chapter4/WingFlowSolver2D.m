
g = Airfoilg();
[p,e,t] = initmesh(g,'hmax',0.5);
A = StiffnessAssembler(p,t,@(x,y) 1);
[R,r] = RobinAssembler2D(p,e,@Kappa1,@gD1,@gN1);
phi = (A + R) \ r;
[phix,phiy] = pdegrad(p,t,phi);
u = -phix';
v = -phiy';

pdegplot(g);
hold on;
% pdeplot(p,e,t,"XYData",phi,"Contour","on","XYStyle","off","FlowData",[u,v]);
pdeplot(p,e,t,"XYData",-sqrt(u .* u + v .* v),"Contour","on");


function z = Kappa1(x,y)
z = 0;
if (x > 29.99)
    z = 1.e+6;
end
end

function z = gD1(x,y)
z = 0;
end

function z = gN1(x,y)
z = 0;
if (x < -29.99)
    z = 1;
end
end