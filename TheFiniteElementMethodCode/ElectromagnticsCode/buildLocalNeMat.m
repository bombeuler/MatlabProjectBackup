function [AK,bK] = buildLocalNeMat(x,y,MU,KAPPA2,fhat)
bK = zeros(3,1);
area = polyarea(x,y);
bb = [y(2) - y(3); y(3) - y(1); y(1) - y(2)]./(2*area);
cc = [x(3) - x(2); x(1) - x(3); x(2) - x(1)]./(2*area);
f = bb*bb' + cc*cc';
m11=(f(3,3)-f(2,3)+f(2,2))/6;
m22=(f(1,1)-f(1,3)+f(3,3))/6;
m33=(f(2,2)-f(1,2)+f(1,1))/6;
m12=(f(3,3)-f(3,1)+2*f(1,2)-f(2,3))/12;
m13=(f(3,2)-2*f(3,1)-f(2,2)+f(2,1))/12;
m23=(f(1,1)-f(1,2)+2*f(3,2)-f(3,1))/12;
MK=[m11 m12 m13; m12 m22 m23; m13 m23 m33].*area;
WK = [1 -1 1;-1 1 -1;1 -1 1]./area;

AK = (1/MU) .* WK + (KAPPA2) .* MK;

bK(1)=dot(fhat,[bb(3); cc(3)]-[bb(2); cc(2)]);
bK(2)=dot(fhat,[bb(3); cc(3)]-[bb(1); cc(1)]);
bK(3)=dot(fhat,[bb(2); cc(2)]-[bb(1); cc(1)]);
bK=bK.*area./3;

end