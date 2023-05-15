function resultZ = problemSolver(p,t)
np = size(p,2);
resultZ = zeros(np,1);
% find nodes lying on edges
edges = findEdgeNodes(p);
usedNodes = setdiff(1:np,edges);
A = StiffnessAssembler(p,t,@(x,y) 1);
b = LoadAssembler2D(p,t,@(x,y) 1);
resultZ(usedNodes) = A(usedNodes,usedNodes)\b(usedNodes);
end