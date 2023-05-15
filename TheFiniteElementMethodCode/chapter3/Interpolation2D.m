function IhF = Interpolation2D()
f = @(x,y) sin(x + y);
hMax = 0.1;

model = createpde(2);
geometryFromEdges(model,@lshapeg);
meshes = generateMesh(model,"hmax",hMax,"GeometricOrder","linear");
p = meshes.Nodes;
t = meshes.Elements;

InF = f(p(1,:),p(2,:));

pdeplot(meshes, ...
        'XYData',InF,'ZData',InF,'Mesh','on');

end