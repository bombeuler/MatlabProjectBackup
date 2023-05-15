function [p,t] = makeMesh(xbl,ybl,width,N)
p = zeros(2,(N+1)*(N+1));
hh = 1/N;
for ii = 0:N
    for jj = 0:N
        index = ii * (N+1) + jj + 1;
        x = xbl + jj*hh*width;
        y = ybl + ii*hh*width;
        p(:,index) = [x,y]';
    end
end

t = zeros(3,2*N*N);

for ii=1:N
    for jj=1:N
        tindex = 2*(jj + (ii-1)*N);
        lb = jj + (ii-1)*(N+1);
        lr = jj + 1 + (ii-1)*(N+1);
        lt = jj + ii*(N+1);
        rt = jj + 1 + ii*(N+1);
        t(:,tindex-1) = [lb lr lt]';
        t(:,tindex) = [lr lt rt]';
    end
end

end