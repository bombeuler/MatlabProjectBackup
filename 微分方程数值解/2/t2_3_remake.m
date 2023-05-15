%% 同时包括内点和外点的
function uSquare = t2_3_remake(N)
% N=2^3;
h=1/N;

M = kron(speye(N+1),gallery('tridiag',2*N+1)) + kron(gallery('tridiag',N+1),speye(2*N+1));

% M = M./(h^2);

pointD = [1:N+1 2*N+1 ((2*N+1)*N +1):(N+1)*(2*N+1) (2*N+2):(2*N+1):(N*(2*N+1)+1) (2*(2*N+1)):(2*N+1):(N*(2*N+1))];
pointG = (N+2):(2*N);

pointIn = setdiff(1:((N+1)*(2*N+1)),pointD);

lengthD = length(pointD);

M(pointD,:) = 0;
M(:,pointD) = 0;
M(pointD,pointD) = speye(lengthD,lengthD);

M(pointG,pointG) = M(pointG,pointG) ./ 2;

% disp(M(pointG,pointG));

% disp(full(M));

F = ones((N+1)*(2*N+1),1) .* (h^2);
F(pointG) = F(pointG) ./ 2;
F(pointD) = 0;

M1 = M;
F1 = ones((N+1)*(2*N+1),1);
F1(pointG) = F1(pointG) ./ 2;
F1(pointD) = 0;

save("M1.mat","M1");
save("F1.mat","F1");

uuM = M\F;

uSquare = reshape(uuM,2*N+1,N+1)';
% disp(uSquare);

% xx = -1:h:1;
% yy = 0:h:1;
% 
% [X,Y] = meshgrid(xx,yy);
% 
% surf(xx,yy,uSquare);

load("uSquare1024.mat","uSquare1024");
disp([errorNorm(uSquare,uSquare1024,h,-1),errorNorm(uSquare,uSquare1024,h,0)]);

end









