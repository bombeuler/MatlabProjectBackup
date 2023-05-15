function uu = deltaProcess(N)
h=2/N;
[p,t] = makeMesh(-1,-1,2,N);
% pdeplot(p,[],[t;ones(1,2*N*N)]);
A = StiffnessAssembler(p,t,1);
np = size(p,2);
midP = (np + 1)/2;

pointB = [1:N+1 ((N+1)*N+1):(N+1)^2 (N+2):(N+1):((N+1)*N+1) (2*N+2):(N+1):(N+1)*N];
pointFree = setdiff(1:np,pointB);

b = sparse(np,1);
b(midP) = 1;
xx = sparse(np,1);
xx(pointFree) = A(pointFree,pointFree)\b(pointFree);

uu = reshape(xx,N+1,N+1)';
% [X,Y] = meshgrid(-1:h:1,-1:h:1);
end
