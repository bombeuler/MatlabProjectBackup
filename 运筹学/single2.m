function [xRes,iterNum] = single2(A,b,c)
M=1e6;
iterNum = 0;
[m,n] = size(A);
basisIndex = 1:m;
x0 = zeros(n,1);
xRes = zeros(n,1);
[basisOrder,A,c] = resortX(A,c);
x0(basisIndex) = b;
sigmaVector = zeros(1,n);
cB = c(basisIndex);
for k = setdiff(1:n,basisIndex)
    sigmaVector(k) = c(k) - cB * A(:,k);
end

while(true)
ADisp = zeros(m,n);
sigmaDisp = zeros(1,n);
ADisp(:,basisOrder) = A;
sigmaDisp(basisOrder) = sigmaVector;
disp([ basisOrder(basisIndex)' b ADisp; -1 -1 sigmaDisp]);
disp(iterNum);
% [~,kk] = max(sigmaVector);

[~,ll] = min(b);
thetaVector = sigmaVector./A(ll,:);
thetaVector(A(ll,:)>=0)=M;
[~,kk] = min(thetaVector);
    
if all(b>=0)
   xRes(basisOrder) = x0; 
   break;
elseif all(A(:,kk)>=0)
   xRes = "unbounded";
   break;
else
    
    alk = A(ll,kk);
    
    
    b(ll) = b(ll)/alk;
    A(ll,:) = A(ll,:)./alk;
    for k = 1:m
        if(k ~=ll)
            b(k) = b(k) -A(k,kk)*b(ll);
            A(k,:) = A(k,:) - A(k,kk).*A(ll,:);
        end
    end
    sigmak = sigmaVector(kk);
    sigmaVector(ll) = -sigmak/alk;
    for k = 1:n
        if(k ~=ll)
            sigmaVector(k) = sigmaVector(k) - sigmak*A(ll,k);
        end
    end
end

basisIndex(ll) = kk;
x0 = zeros(n,1);
x0(basisIndex) = b;

iterNum = iterNum + 1;
end
end

function [basisOrder,ASort,cSort] = resortX(A,c)
[m,n] = size(A);
v=nchoosek(1:n,m);
for i=1:size(v,1)
    if A(:,v(i,:))==eye(m)
        basisOrder=[v(i,:) setdiff(1:n,v(i,:))];
    end
end
cSort = c(basisOrder);
ASort = A(:,basisOrder);

end
