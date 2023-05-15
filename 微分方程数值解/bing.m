% Define geometry
g = [1;0;0;1]; % Circle with radius 1 and center (0,0)
gd = [g [1;0.5;0.5;0.5]]; % Circle with radius 0.5 and center (0.5,0.5)
sf = 'C1-C2'; % Set formula for geometry
ns = char('C1','C2')'; % Assign names to geometry
disp(decsg(gd,sf,ns));
[p,e,t] = initmesh(decsg(gd,sf,ns)); % Create geometry and mesh

% Define boundary conditions
b = pderect([-1 1 -1 1],'C1'); % Outer boundary of C1
b(:,5) = 0; % Set all outer boundary segments to Neumann condition
b(3,3) = 100; % Set Dirichlet condition of 100 V on segment 3
b(3,4) = -100; % Set Dirichlet condition of -100 V on segment 4
c = pdeintrp(p,t,pdecirc(0.5,0.5,0.5)); % Inner boundary of C2
c(:,5) = 0; % Set all inner boundary segments to Neumann condition

% Solve PDE
u = assempde(b,p,e,t,c,1,0,0); % Assemble and solve PDE
pdeplot(p,e,t,'xydata',u,'zdata',u,'mesh','on') % Plot solution and mesh
pdetool % Open graphical interface