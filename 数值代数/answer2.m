format short;

%% 定义A,b,L,U,D,以及精确解x
A0=[16,2,3,1;1,10,4,3;3,1,15,2;1,2,4,18];
b0=[15;1;-40;61];
A=A0'*A0;
b=A0'*b0;
disp(A);
disp(b);
L=zeros(4,4);
U=zeros(4,4);
for i=1:4
    for j=1:4
        if i>j
            L(i,j)=A(i,j);
            U(j,i)=A(j,i);
        end
    end
end
D=A-L-U;
x=[3887/2948;568/4085;-2337/671;1813/445];

%% Jacobi迭代法
x1=zeros(4,1);
x2=-inv(D)*(L+U)*x1+D\b;
k=1;
p1=vrho(-inv(D)*(L+U));
disp(p1);
while(norm(x1-x2,inf)>=10^(-5))
    x1=x2;
    x2=-inv(D)*(L+U)*x2+D\b;
    k=k+1;
end
disp(k);
disp(x2);
disp(norm(x-x2));

%% Gauss-Seidel迭代法
x1=zeros(4,1);
x2=-inv(D+L)*U*x1+(D+L)\b;
k=1;
p2=vrho(-inv(D+L)*U);
disp(p2);
while(norm(x1-x2,inf)>=10^(-5))
    x1=x2;
    x2=-inv(D+L)*U*x2+inv(D+L)*b;
    k=k+1;
end
disp(k);
disp(x2);
disp(norm(x-x2));

%% JOR迭代法
lamda=eig(-inv(D)*(L+U));
w=2/(2-max(lamda)-min(lamda));
disp(w);
p3=vrho(eye(4)-w*inv(D)*A);
disp(p3);
x1=zeros(4,1);
x2=(eye(4)-w*inv(D)*A)*x1+w*inv(D)*b;
k=1;
while(norm(x1-x2,inf)>=10^(-5))
    x1=x2;
    x2=(eye(4)-w*inv(D)*A)*x2+w*inv(D)*b;
    k=k+1;
end
disp(k);
disp(x2);
disp(norm(x-x2));f

%%SOR迭代法
w=2/(1+sqrt(1-(vrho(-inv(D)*(L+U)))^2));
disp(w);
p4=vrho(inv(D+w*L)*((1-w)*D-w*U));
disp(p4);
x1=zeros(4,1);
 x2=inv(D+w*L)*(((1-w)*D-w*U)*x1+w*b);
 k=1;
while(norm(x1-x2,inf)>=10^(-5))
    x1=x2;
    x2=inv(D+w*L)*(((1-w)*D-w*U)*x2+w*b);
    k=k+1;
end
disp(k);
disp(x2);
disp(norm(x-x2));

function p=vrho(A)
p=max(abs(eig(A)));
end