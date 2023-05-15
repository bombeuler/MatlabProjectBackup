% h=0.01;
function uu=t1_1(h,flag)
N=1/h;
tt=0:h:1;
uu=zeros(size(tt));
U=exp(-2.*tt);
uu(1)=1;
for i=1:N
    if flag==1
        uu(i+1)=(1-2*h)*uu(i);
    elseif flag==2
        uu(i+1)=uu(i)/(1+2*h);
    else
        if i>1
            uu(i+1)=uu(i-1)-4*h*uu(i);
        end
    end
end
% plot(tt,uu);
plot(tt,abs(uu-U));
end