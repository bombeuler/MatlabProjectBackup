function PhF = L2Projector2D()
f = @(x,y) (sin(x + y))  + exp(x +y);
hMax = 0.1;

model = createpde(2);
geometryFromEdges(model,@lshapeg);
meshes = generateMesh(model,"hmax",hMax,"GeometricOrder","linear");
p = meshes.Nodes;
t = meshes.Elements;
M = MassAssembler2D(p, t);
b = LoadAssembler2D(p, t, f);
PhF = M \ b;
spy(M);
disp(size((PhF)));
% pdesurf(p, t, Pf);
% pdeplot(meshes, ...
%         'XYData',PhF,'ZData',PhF,...
%         'FaceAlpha', 0.4, 'ColorMap','jet','Mesh','on');
end