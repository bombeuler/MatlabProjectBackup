% h=0.01;
function uu=t1_2(h,flag)
N=1/h;
tt=0:h:1;
uu=zeros(1,N+1);
U=1./(1-tt);
uu(1)=1;
for i=1:N
    if flag==1
        uu(i+1)=uu(i)+h*(uu(i))^2;
    elseif flag==2
%         res=solve(x-uu(i)+h.*x.^2,x);
%         uu(i+1)=res(res>0);
        uu(i+1)=(1+sqrt(1-4*h*uu(i)))/(2*h);
    else
        if i>1
            uu(i+1)=uu(i-1)+2*h*(uu(i))^2;
        end
    end
end
% plot(tt,uu);
plot(tt,abs(uu-U));
end