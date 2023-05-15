%|********************************************************************%
%| This demo is to use the central finite difference scheme to solve  %  
%| the wave equation with zero Dirichlet boundary condition           %
%|                                                                    %  
%|  u_tt = u_xx         for 0<x<1                                     %    
%|  u(0,t) = u(1,t) = 0 on the boundary                               %
%|  u(x,0) = sin(pi*x), u_t(x,0) = sin(pi*x), at t = 0                %               
%|                                                                    %
%|  The exact solution is u(x,t) = sin(pi*x)(cos(pi*t)+sin(pi*t)/pi)  %
%|********************************************************************%
%| This short matlab code is as demo code in my class 2022 Spring.    %
%| If you find mistakes or have a better idea in the implementation,  %
%| please send an email to Huadong GAO (Email:huadong@hust.edu.cn)    %
%*********************************************************************%

clear;

T = 1.5;

for kkk =1:4
    
    N = 4*2^kkk;    % spatial grid number
    h =   1.0/N;    % mesh size
    x =   0:h:1;    % grid
    dt =      h;    % time step

    % generate [-1 2 -1] matrix. We only compute interior grids, as B.C. is 0.
    A = gallery('tridiag',N-1); 
    
    uh0= u_initial(x);
    uh1= uh0 + dt*ut_initial(x) + dt*dt/2*(-1/h/h)*[0 (A*uh0(2:(end-1))')' 0];
    uh2= zeros(size(uh0));
    
    r = dt/h; % CFL number
    
    tc = 2*dt;
    while (tc < T+1e-13)

        % the main computation
        uh2(2:(end-1)) = -uh0(2:(end-1)) + 2*uh1(2:(end-1)) - r*r*(A*uh1(2:(end-1))')';
        
        if(kkk==4) % if the finest mesh is used, then we plot the wave
            clf
            subplot(1,3,1)
            plot(x,u_initial(x),'k-'), hold on
            plot(x,uh2,'b-o','linewidth',2),
            title(['numerical solution at t = ' num2str(tc)])
            hold off
            axis([-0.1 1.1 -1.1 1.1])
            subplot(1,3,2)
            plot(x,u_initial(x),'k-'), hold on
            plot(x,u_ex(x,tc),'r-'),
            title(['exact solution at t = ' num2str(tc)])
            hold off
            axis([-0.1 1.1 -1.1 1.1])
            subplot(1,3,3)
            plot(x,uh2-u_ex(x,tc),'k-d'),
            title(['error at t = ' num2str(tc)])
            hold off
            drawnow
        end

        % update some data
        tc         = tc + dt;
        uh0        = uh1;
        uh1        = uh2;
        
    end

    % compute errors
    tc = tc -dt;    
    error = norm(u_ex(x,tc)-uh2)*sqrt(h);
    disp(num2str([N tc error]))

end


function y = u_initial(x) % initial solution u
    y = sin(pi*x);
end

function y = ut_initial(x) % initial solution u_t
    y = sin(pi*x);
end

function y = u_ex(x,t) % exact solution
    y = sin(pi*x).*(cos(pi*t)+sin(pi*t)/pi);
end
