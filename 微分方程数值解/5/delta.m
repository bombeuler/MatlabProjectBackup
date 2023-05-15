clc,clear

g=Rectg(0,0,1,1);
%��ȷ��
N=512;
h=1/N;
[p,e,t]=poimesh(g,N,N);
np=size(p,2);
u=zeros(np,1);
np = size(p,2);
nt = size(t,2);
A=assema(p,t,1,0,0);
point=N^2/2+N+1;
b=sparse(point,1,1,np,1);
fixed =unique([e(1,:) e(2,:)]);
free =setdiff (1: np, fixed);
u(free)=A(free,free)\b(free);
realu=u;



%��������
N=512;
k=log2(N)-2;
Lerror=zeros(k,1);
maxerror=zeros(k,1);
Lratio=zeros(k-1,1);
maxratio=zeros(k-1,1);

for j=1:k
    N=N/2;
    h=2/N;
    [p,e,t]=poimesh(g,N,N);
    figure(j);
    np=size(p,2);
    nt=size(t,2);
    u=zeros(np,1);
    A=assema(p,t,1,0,0);
    
    point=N^2/2+N+1;
    b=sparse(point,1,1,np,1);
    fixed=unique([e(1,:) e(2,:)]);
    free=setdiff (1: np, fixed);
    u(free)=A(free,free)\b(free);
    [X,Y ] = meshgrid(-1:h:1,-1:h:1);
    mesh(X,Y,reshape(u,N+1,N+1)');
    
    %��ϸ�����ϵľ�ȷ��ӳ�䵽��������
    number=1:2:2*N+1;
    row=1:(N+1)^2;
    col=number;
    for i=1:N
        number=number+4*N+2;
        col=[col,number];
    end
    realu=realu(col);
    
    
    
    % ���
    error=u-realu;
    Lerror(j)=h*norm(abs(error),2);
    maxerror(j)=max(abs(error));
end

%�����ı�ֵ  ��������
for i=1:k-1
    maxratio(i)=maxerror(i+1)/maxerror(i);
    Lratio(i)=Lerror(i+1)/Lerror(i);
end
maxratio=log2(maxratio);
Lratio=log2(Lratio);






