function b = LoadAssembler2D(p, t, f)
nP = size(p, 2);
nT = size(t, 2);
b = zeros(nP, 1);
for k = 1:nT
    pointK = t(1:3, k);
    xx = p(1, pointK);
    yy = p(2, pointK);
    areaK = polyarea(xx, yy);
    bK = [f(xx(1), yy(1));
          f(xx(2), yy(2));
          f(xx(3), yy(3))] .* areaK ./ 3;
    b(pointK) = b(pointK) + bK;
end
end