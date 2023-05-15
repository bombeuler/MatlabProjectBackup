% clear;clc;
g = [2     2     2     2;
    -1     1     1    -1;
     1     1    -1    -1;
     1     1     0     0;
     1     0     0     1;
     0     0     0     0;
     1     1     1     1];
% model = createpde(2);
% gg = geometryFromEdges(model,g);
% meshes = generateMesh(model,"hmax",0.2,"GeometricOrder","linear");
% [p,e,t] = meshToPet(meshes);

kmax = 9;
errors = zeros(1,kmax);

for kk = 1:kmax

[p,e,t] = poimesh(g,2^kk,2^(kk-1));

resultZ = problemSolver(p,t);

% for ii=1:1
% 
% eta = pdejmps(p,t,1,0,1,resultZ,1,1,1);
% tol = 0.8*max(eta);
% elements = find(eta > tol);
% 
% [p,e,t] = refinemesh(g,p,e,t,elements',"regular");
% resultZ = problemSolver(p,t);
% end

eta2 = pdejmps(p,t,1,0,1,resultZ,1,1,1);
errors(kk) = (sum(eta2.^2));

end

plot(1:kmax,errors);
% figure(1);
% pdeplot(p,e,t);
% title(strcat("nt=",num2str(size(t,2))));
% figure(2);
% pdeplot(p,e,t,'XYData',resultZ,'ZData',resultZ);