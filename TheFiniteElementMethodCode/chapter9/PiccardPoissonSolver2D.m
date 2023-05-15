g = [2     2     2     2
     0     1     1     0
     1     1     0     0
     0     0     1     1
     0     1     1     0
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'hmax',0.05); % create mesh
xi=zeros(size(p,2),1); % initial zero guess

Afcn = @(u) 0.001+ u .^2;
Ffcn = @(x,y) x .^ 0;

for k=1:1 % non-linear loop
   [A,bb] = PiccardStiffnessLoadAssembler(p,t,xi,Afcn,Ffcn);
    d=A\bb; % solve for correction
    xi=xi+d; % update solution
    fprintf('epoch %d : |d|=%.18f,|r|=%.18f\n', k,norm(d), norm(bb));
end

pdesurf(p,t,xi);