
%**********************************************************************%
%| This demo is to use the 5 points central finite difference scheme   %
%| for the Poisson equation with mixed boundary condition              %
%|               - \Delta u=f(x,y),  x \in \Omega=[-1,1]x[0,1]         %
%| (\partial u)/(partial n)=0,       x on right half of bottom boundary%
%|                        u=0 on other parts of the boundary.          %
%|                   f(x,y)=1 is used in the test                      %
%|*********************************************************************%
%| This short matlab code is as demo code in my class 2022 Spring.     %
%| If you find mistakes or have a better idea in the implimentation,   %
%| please send an email to Huadong GAO (Email:huadong@hust.edu.cn)     %
%**********************************************************************%

for kkk = 0

    h  =            1/8/2^kkk; % a uniform mesh size is used, check convergence rates by refining
    xh =               -1:h:1; % nodes in x-direction
    yh =                0:h:1; % nodes in y-direction
    Nx =           length(xh); % node number in x-direction
    Ny =           length(yh); % node number in y-direction
    A  = kron(speye(Ny),gallery('tridiag',Nx)) + kron(gallery('tridiag',Ny),speye(Nx)); % generate matrix
    A  =              1/h/h*A; % % generate matrix
    
    itot    =                 Nx*Ny; % total  nodes index
    ibd_C   = [1,Nx,itot-Nx+1,itot]; % corner nodes index
    ibd_D   = [2:((Nx+1)/2) (Nx+1):Nx:(itot-2*Nx+1) (itot-Nx+2):(itot-1) (2*Nx):Nx:(itot-Nx)]; % Dirichlet nodes
    ibd_all = [2:(Nx-1)     (Nx+1):Nx:(itot-2*Nx+1) (itot-Nx+2):(itot-1) (2*Nx):Nx:(itot-Nx)]; % all nodes on boundary
    ibd_N   = [((Nx+1)/2+1):(Nx-1)];                                                           % Neumann nodes
    
    % generate right-hand-side
    [xxh,yyh]  =    meshgrid(xh,yh); %  position 
    F          =     f_rhs(xxh,yyh); %  take values
    F          = reshape(F',itot,1); % 
    
    % modify the right half bottom boundary with Neumann B.C. 
    for k = 1:length(ibd_N)
        A(ibd_N(k),(ibd_N(k)-1):(ibd_N(k)+1)) = A(ibd_N(k),(ibd_N(k)-1):(ibd_N(k)+1))/2;
    end
    
    % modify 4 corner nodes with  Dirichlet B.C.
    A(ibd_C,:) = A(ibd_C,:)*0;
    A(:,ibd_C) = A(:,ibd_C)*0;
    for k=1:4
        A(ibd_C(k),ibd_C(k)) = 1;
        F(ibd_C(k)         ) = 0;
    end
    
    % modify 4 boundary nodes with  Dirichlet B.C.
    A(ibd_D,:) = A(ibd_D,:)*0;
    A(:,ibd_D) = A(:,ibd_D)*0;
    for k = 1:length(ibd_D)
        A(ibd_D(k),ibd_D(k)) = 1;
        F(ibd_D(k)         ) = 0;
    end
    
    uh  = A\F; % solve linear system
    
    % weigth of each node in discrete form integral 
    b_vec          = ones(Nx*Ny,1);
    b_vec(ibd_C)   =           1/4;
    b_vec(ibd_all) =           1/2;
    
    load('uh1024.mat')
    ntab = h*1024;
    ue   = uh1024(1:ntab:end,1:ntab:end); % exact solution is computed with
    ue   =          reshape(ue', itot,1); % a 2048x1024 uniform mesh
    errl2   =        sqrt(h*h*b_vec'*(uh - ue).^2); % l2-error
    errlinf =                     norm(ue-uh, inf); % l-inf error
    
    disp('   l2 error:                      l-inf error : ')
    format longe
    disp([errl2, errlinf])

end


uh = reshape(uh(1:itot),Nx,Ny)';
colormap('jet')
subplot(1,2,1)
surf(xxh, yyh, uh)
shading interp 
%axis([-1.1 1.1 -0.1 1.1 -0.1 0.2])
title('Numerical solution of 5-point central FDM')

subplot(1,2,2)
contourf(xxh,yyh,uh,64)
title('Numerical solution of 5-point central FDM')

% 定义 right-hand-side 右端函数
function y = f_rhs(x1,x2)
    y = 1 + 0*x1;
end