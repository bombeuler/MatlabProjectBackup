% 求离散范数的函数
%0为0，c为2
function res = disperse_norm2(u,flag)
Nlist=size(u)-[1,1];
h=1./Nlist;
hs=prod(h);
rho_m=zeros(Nlist);

for i=1:Nlist(1)+1
    for j=1:Nlist(2)+1
        rho_m(i,j)=rho(i,j);
    end
end

switch flag
    case 2
%         res=max(abs(u),[],'all');
           res=max(max(abs(u)));
    case 0
%         res=sqrt(sum(u.^2.*rho_m,'all'));
           res=sqrt(sum(sum(u.^2.*rho_m)));
end

    function ss=rho(i,j)
        if (i==1 && j==1) || (i==1 && j==Nlist(2)+1) || ( i==Nlist(1)+1 &&j==1) || (i==Nlist(1)+1 &&j==Nlist(2)+1)
            ss=hs/4;
        elseif i==1 || i==Nlist(1)+1 || j==1 || j==Nlist(2)+1
            ss=hs/2;
        else
            ss=hs;
        end
    end

end