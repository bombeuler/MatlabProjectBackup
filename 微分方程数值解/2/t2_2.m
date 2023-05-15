
h=0.1;
N=2/h;
Nm=1/h;

xx1=-1:h:1;
yy1=-1:h:0;

xx2=0:h:1;
yy2=0:h:1;

points_b=(N+1)*Nm+Nm;
points_t=(Nm+1)^2;
points=points_b+points_t;

Mb=zeros(points_b);
Mt=zeros(points_t);

A1=gallery('poisson',N-1);
A2=gallery('poisson',Nm-1);
O1=zeros(N-1,1);
O2=zeros(Nm-1,1);

%处理Mt 

Mt(1:Nm+1,1:Nm+1)=eye(Nm+1);
Mt(Nm+2,Nm+2)=1;
Mt(2*(Nm+1),2*(Nm+1))=1;
Mt(Nm+3:2*Nm+1,1:2*(Nm+1))=[O2 A2(1:Nm-1,1:Nm-1) O2 O2 A2(1:Nm-1,Nm:2*(Nm-1)) O2];

for i=3:Nm-2
    pre_i=(i-1)*(Nm+1)+1;
    next_i=i*(Nm+1);
    pre_a=(i-3)*(Nm-1);
    mid_a=(i-2)*(Nm-1);
    next_a=(i-1)*(Nm-1);
    i_a=i*(Nm-1);
    pm_a=pre_a+1:mid_a;
    mn_a=mid_a+1:next_a;
    ni_a=next_a+1:i_a;
    
    %赋值
    Mt(pre_i,pre_i)=1;
    Mt(next_i,next_i)=1;
    Mt(pre_i+1:next_i-1,(i-2)*(Nm+1)+1:(i+1)*(Nm+1))=[O2 A2(mn_a,pm_a) O2 O2 A2(mn_a,mn_a) O2 O2 A2(mn_a,ni_a) O2];
end

Mt(points_t,Nm*(Nm+1)+1)=1;
Mt(points_t,points_t)=1;

Mt(Nm*(Nm+1)+2:Nm*(Nm+1)+Nm,(Nm-1)*(Nm+1)+1:(Nm+1)^2)=[O2 A2((Nm-3)*(Nm-1)+1:(Nm-2)*(Nm-1),(Nm-3)*(Nm-1)+1:(Nm-2)*(Nm-1)) O2 O2 A2((Nm-2)*(Nm-1)+1:(Nm-1)^2,(Nm-2)*(Nm-1)+1:(Nm-1)^2) O2];
Mt(Nm*(Nm+1)+2:Nm*(Nm+1)+Nm,Nm*(Nm+1)+2:Nm*(Nm+1)+Nm)=-eye(Nm-1);

% 处理Mb 

% 前Nm行

Mb(1:Nm,1:Nm)=eye(Nm);

%Nm+1 ~ 2 Nm+1 ~Nm+N+1
Mb(Nm+1,Nm+1)=1;
Mb(Nm+N+1,Nm+N+1)=1;
Mb(Nm+2:Nm+N,Nm+1:Nm+2*N+2)=[O1 A1(1:N-1,1:N-1) O1 O1 A1(1:N-1,N:2*N-2) O1];

for i=2:Nm-2
    pre_i=i*(N+1)+1+Nm;
    next_i=(i+1)*(N+1)+Nm;
    pre_a=(i-1)*(N-1);
    mid_a=i*(N-1);
    next_a=(i+1)*(N-1);
    i_a=(i+2)*(N-1);
    pm_a=pre_a+1:mid_a;
    mn_a=mid_a+1:next_a;
    ni_a=next_a+1:i_a;
    
    %赋值
    Mb(pre_i,pre_i)=1;
    Mb(next_i,next_i)=1;
    Mb(pre_i+1:next_i-1,(i-2)*(N+1)+1+Nm:(i+1)*(N+1)+Nm)=[O1 A1(pm_a,pm_a) O1 O1 A1(pm_a,mn_a) O1 O1 A1(pm_a,ni_a) O1];
end


Mb(Nm+(Nm-2)*(N+1)+1,Nm+(Nm-2)*(N+1)+1)=1;
Mb(Nm+(Nm-1)*(N+1),Nm+(Nm-1)*(N+1))=1;
Mb(Nm+(Nm-2)*(N+1)+2:Nm+(Nm-2)*(N+1)+N,(Nm-2)*(N+1)+Nm+1:Nm*(N+1)+Nm)=[O1 A1((Nm-2)*(N-1)+1:(Nm-1)*(N-1),(Nm-2)*(N-1)+1:(Nm-1)*(N-1)) O1 O1 A1((Nm-1)*(N-1)+1:Nm*(N-1),(Nm-1)*(N-1)+1:Nm*(N-1)) O1];

Mb(Nm+(Nm-1)*(N+1)+1:Nm+Nm*(N+1),Nm+(Nm-1)*(N+1)+1:Nm+Nm*(N+1))=eye(N+1);

M=blkdiag(Mt,Mb);

Ft=ones(Nm+1);
Ft(:,1)=0;
Ft(1,2:Nm+1)=0;
Ft(Nm+1:2:Nm+1)=0;

Fm=zeros(Nm,1);

Fb=ones(N+1,Nm);
Fb(1,:)=0;
Fb(N+1,:)=0;
Fb(:,Nm)=0;

Fl=[reshape(Ft,(Nm+1)^2,1);Fm;reshape(Fb,(N+1)*Nm,1)];
disp([size(Mt),size(Mb)]);
disp([size(Fl);size(M)]);

ul=M\Fl;
disp(ul')





