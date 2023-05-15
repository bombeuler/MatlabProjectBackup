function M = MassAssembler2D(p,t)
nP = size(p,2);
nT = size(t,2);
M = sparse(nP,nP);
for k = 1:nT
    pointK = t(1:3,k);
    areaK = polyarea(p(1,pointK), p(2,pointK));
    MK = [2 1 1;
          1 2 1;
          1 1 2] .* areaK ./ 12;
    M(pointK,pointK) = M(pointK,pointK) + MK;

end
end