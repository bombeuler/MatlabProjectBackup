function [A,bb] = PiccardStiffnessLoadAssembler(p,t,u,Afcn,Ffcn)
np = size(p,2);
nt = size(t,2);
A = sparse(np,np);
for K = 1:nt
    loc2glb = t(1:3,K);
    x = p(1,loc2glb);
    y = p(2,loc2glb);
    [area, b, c ] = HatGradients(x,y);
    amean = mean(Afcn(u(loc2glb)));
    disp(amean);
    AK = amean * (b * b' + c * c') * area;
    A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK;
end
bb = LoadAssembler2D(p,t,Ffcn);
end

