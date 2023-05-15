% [p,e,t] = initmesh('squareg','hma',0.05); % mesh
% x = p(1,:); y = p(2,:); % node coordinates
% pif = x'*y;% nodal values of interpolant
% pdesurf(p,t,pif) % plot interpolant
xx=0:0.02:1;
plot([0,1/3,2/3,1],[2,28/9,31/9,3],xx,(4 + 9.*xx - 7.*xx.^2)./2);