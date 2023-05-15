%|********************************************************************%
%|         This demo is to use the leap-frog scheme to solve          %
%|     the advection equation with periodic  boundary condition       %
%|                                                                    %  
%|  u_t + u_x = 0       for 0<x<1                                     %
%|                                                                    %
%|  The intial solution is a smooth hat function                      %
%|        u0(x,0) = exp(-60*pi*(x-0.3).^2)                            %
%|                                                                    %
%|  The exact solution is u(x,t) = u0(x-t,0)                          %
%|********************************************************************%
%| This short matlab code is as demo code in my class 2022 Spring.    %
%| If you find mistakes or have a better idea in the implementation,  %
%| please send an email to Huadong GAO (Email:huadong@hust.edu.cn)    %
%|********************************************************************%

T  =   1.5;
N  =   128;
h  = 1.0/N;
x  = 0:h:1;
dt =     h;
r  =  dt/h;

uh0 =    u_initial(x); %  Inital time step
uh1 = u_initial(x-dt); %  We take the exact solution at the first step
uh2 = zeros(size(uh0));

tc = 2*dt;
while (tc < T+1e-13)

    uh2(2:end) = uh0(2:end) - r*([uh1(3:end) uh1(1)]-uh1(1:(end-1)));
    uh2(1)     = uh2(end);

    clf    
    plot(x,u_initial(x),'r-*'), hold on    
    hold on, plot(x,uh2,'b-o','linewidth',2), hold off
    axis([-0.1 1.1 -0.1 1.1])
    drawnow

    tc         = tc + dt;
    uh0        = uh1;
    uh1        = uh2;
    
end

tc = tc -dt;
disp(tc)

function y = u_initial(x)
    y = exp(-60*pi*(x-0.3).^2);
end