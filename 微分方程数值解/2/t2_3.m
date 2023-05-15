%% 只包括内点的
function uu=t2_3(N)
% N=2^10;
h=1/N;

M_in=gallery('poisson',2*N-1);

%% M
M=sparse(2*N*(N-1),2*N*(N-1));

% 1
M(1:N-1,1:N-1)=gallery('tridiag',N-1,-1/2,2,-1/2);
M(1:N-1,2*N:3*N-2)=-speye(N-1);

% 2
M(2*N:3*N-2,1:N-1)=-speye(N-1);
M(N:2*N*(N-1),N:2*N*(N-1))=M_in(1:(2*N-1)*(N-1),1:(2*N-1)*(N-1));

% M=M/(h^2);

%% F
% F1=ones(1,N-1);
% F2=ones(N-1,2*N-1);
Fl=h^2.*ones(2*N*(N-1),1);
Fl(1:N-1) = Fl(1:N-1)/2;

% spy(M-M');

M2 = M;
F2=ones(2*N*(N-1),1);
F2(1:N-1) = F2(1:N-1)/2;
save("M2.mat","M2");
save("F2.mat","F2");

ul=M\Fl;

% disp(M*ul .*(N^2));

%% uu
uu=zeros(N+1,2*N+1);
uu(1,N+2:2*N)=ul(1:N-1)';
uu(2:N,2:2*N)=reshape(ul(N:end),2*N-1,N-1)';

% disp(uu);

load("uTrue.mat","uTrue");

errorC = errorNorm(uu,uTrue,h,-1);
error0 = errorNorm(uu,uTrue,h,0);
disp([errorC,error0]);

% % plot
% xx=-1:h:1;
% yy=0:h:1;
% [X,Y]=meshgrid(xx,yy);
% % disp([X;Y]);
% 
% figure(1);
% mesh(X,Y,uu);

end




