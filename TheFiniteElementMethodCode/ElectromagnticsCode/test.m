dl = decsg([2 4 -5 5 5 -5 -5 -5 5 5]');
[p,e,t] = poimesh(dl,2,2);
pdeplot(p,e,t,Nodelabels="on");
[ee,te,t2e] = buildEdges(p,t(1:3,:));
OMEGA = 2*pi;
EPSILON = 1;
MU = 1;
SIGMA = 0;

KAPPA2 = 1i*OMEGA .* SIGMA - OMEGA^2 .* EPSILON; 

np = size(p,2);
ne = size(ee,2);
nt = size(te,2);

A = sparse(ne,ne);
b = sparse(ne,1);

jN = exp(1i*OMEGA/6);

fixed = [];
neuP = [];

for ek = 1:ne
    meanxy = mean(p(:,ee(:,ek)),2);
    xmid = meanxy(1);
    ymid = meanxy(2);

    if(xmid == -5 || xmid == 5)
        fixed = [fixed ek];
    elseif (ymid == -5)
        neuP = [neuP ek];
    end
end

for kk = 1:nt
    x = p(1,te(1:3,kk));
    y = p(2,te(1:3,kk));
    [AK,bK] = buildLocalNeMat(x,y,MU,KAPPA2,[0;0]);
    A(te(4:6,kk),te(4:6,kk)) = A(te(4:6,kk),te(4:6,kk)) + AK;
    b(te(4:6,kk)) = b(te(4:6,kk)) + bK;
end


for ek = 1:neuP
    tk =max(t2e(:,ek));
    meanxy = mean(p(:,ee(:,ek)),2);
    x = p(1,te(1:3,ek));
    y = p(2,te(1:3,ek));
    area = polyarea(x,y);

    if (ek == te(4,ek))
        b(ek) = b(ek)+1i*OMEGA*MU*jN*(y(1)-meanxy(2)+meanxy(1)-x(1))/2/area;
    elseif (ek == te(5,ek))
        b(ek) = b(ek)+1i*OMEGA*MU*jN*(meanxy(2)-y(2)+x(2)-meanxy(1))/2/area;
    else
        b(ek) = b(ek)+1i*OMEGA*MU*jN*(y(3)-meanxy(2)+meanxy(1)-x(3))/2/area;
    end

end

free = setdiff(1:ne,fixed);


