%B=[3 2 -1 2;1 4 0 2;2 1 2 -1;1 1 -1 3];
%c=[-2;1;3;4];
%y=chase([-1;-1;-1;-1],[2;2;2;2;2],[-1;-1;-1;-1],[1;0;0;0;0]);
%[y,L,U]=lux(B,c);
%disp(y);
%disp(Q);
%disp(R);
% A=zeros(6,4);
% x=[0.1,0.3,0.5,0.6,0.7,0.9];
% b=[0.61;0.92;1.12;1.52;1.47;2.04];
%  for i=1:6
%     A(i,1)=1;
%     A(i,2)=x(i);
%     A(i,3)=sin(x(i));
%     A(i,4)=exp(x(i));
%  end
   % disp(A);
%     a=(A' * A)\(A' *b);
%     disp(a);
format short;
% A=[10,7,8,7;7,5,6,5;8,6,10,9;7,5,9,10];
% A=A+0.01.*A;
% b=[32;23;33;31];
% b=b+[0.01;-0.01;0.01;-0.01];
% [x,Q,R]=lux(A,b);
% e=norm(x-[1;1;1;1],inf)/norm([1;1;1;1],inf);
% disp(e);
% disp(A);
% disp(x);
% disp(Q);
% disp(R);
x=[0.1;0.2;0.3;0.4;0.5];
A=[[1;1;1;1;1],x];
disp(A);
b=[log(80.4);log(53.9);log(36.1);log(24.2);log(16.2)];
disp(b);
a=(A'*A)\(A'*b);
disp(a);
% [y,Q,R]=zqrfact(A,b);
% disp(y);
% disp(Q);
% disp(R);


