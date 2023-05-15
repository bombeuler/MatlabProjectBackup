g = [2     2     2     2
     0     1     1     0
     1     1     0     0
     0     0     1     1
     0     1     1     0
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'hmax',0.05); % create mesh
xi=zeros(size(p,2),1); % initial zero guess

for k=1:9 % non-linear loop
    [J,r]=JacResAssembler2D(p,e,t,xi,@Afcn,@Ffcn);
    d=J\r; % solve for correction
    xi=xi+d; % update solution
    fprintf('epoch %d : |d|=%.18f,|r|=%.18f\n', k,norm(d), norm(r));
end

pdesurf(p,t,xi);

function z = Afcn(u)
z=0.125+u.^2;
end

function z = Ffcn(x,y)
z=x.^0; % =1
end