function [AK,bK] = RotRot(x,y,sigmas,mu,kappa,fhat);
[area,b,c]=HatGradients(x,y);
f=b*b'+c*c';
len=zeros(3,1); % edge lengths
len(1)=sqrt((x(3)-x(2))^2+(y(3)-y(2))^2);
len(2)=sqrt((x(1)-x(3))^2+(y(1)-y(3))^2);
len(3)=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
len=len.*sigmas; % edge lengths times signs
m11=(f(3,3)-f(2,3)+f(2,2))/6;
m22=(f(1,1)-f(1,3)+f(3,3))/6;
m33=(f(2,2)-f(1,2)+f(1,1))/6;
m12=(f(3,1)-f(3,3)-2*f(2,1)+f(2,3))/12;
m13=(f(3,2)-2*f(3,1)-f(2,2)+f(2,1))/12;
m23=(f(1,2)-f(1,1)-2*f(3,2)+f(3,1))/12;
MK=[m11 m12 m13; m12 m22 m23; m13 m23 m33]*area;
WK=ones(3)/area;
AK=(WK/mu+kappa^2*MK).*(len*len');
bK=zeros(3,1);
bK(1)=dot(fhat,[b(3); c(3)]-[b(2); c(2)]);
bK(2)=dot(fhat,[b(1); c(1)]-[b(3); c(3)]);
bK(3)=dot(fhat,[b(2); c(2)]-[b(1); c(1)]);
bK=bK.*len*area/3;