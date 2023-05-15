% v=[1 -2 1];
% m=[2 3];
% A=diag(v)
% B=diag(v,1)
% C=diag(v,-1)
% D=diag(m,2)
A1=ones(2);
A2=ones(3);
A=blkdiag(A1,A2);
disp(A);
