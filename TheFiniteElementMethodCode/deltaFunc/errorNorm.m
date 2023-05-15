function errorResult = errorNorm(uu,uTrue,N,flag,fixed)
NTrue = 2^9;
h = 2/N;
transformRatio = NTrue/N;
xx = zeros(N+1,N+1);
tao = ones(N+1,N+1);

tao(1,1) = 1/4;
tao(1,N+1) = 1/4;
tao(N+1,1) = 1/4;
tao(N+1,N+1) = 1/4;

tao(1,2:N) = 1/2;
tao(N+1,2:N) = 1/2;
tao(2:N,1) = 1/2;
tao(2:N,N+1) = 1/2;

if (nargin < 5)

for ii = 1:N+1
    for jj = 1:N+1
        xx(ii,jj) = uu(ii,jj) - uTrue((ii - 1)*transformRatio + 1,(jj - 1)*transformRatio + 1);
%         disp([ii,jj;(ii - 1)*transformRatio + 1,(jj - 1)*transformRatio + 1]);
    end
end

else
 xx = uu - uTrue;
end



if flag == -1
    errorResult = max(abs(xx),[],'all');
elseif flag == 0
    errorResult = sqrt(sum((h^2) .* tao .* (abs(xx) .^2),'all'));
end
end