function problem6(xx, f)
M = assembleMass1D(xx);
b = assembleLoad1D(xx, f);
a = M \ b;
plot(xx, [0 a' 0]);
hold on;
plot(xx, xx .^ 2 - xx, '-r');
end

function M = assembleMass1D(xx)
n = length(xx) - 2;
M = zeros(n,n);

for k = 2:n
    h = xx(k+1) - xx(k);
    M(k-1, k-1) = M(k-1, k-1) + 1 / h + h / 3;
    M(k-1, k) = M(k-1, k) + h / 6 - 1 / h;
    M(k, k-1) = M(k, k-1) + h / 6 - 1 / h;
    M(k, k) = M(k, k) + 1 / h + h / 3;
end
h1 = xx(2) - xx(1);
hEnd = xx(n+2) - xx(n+1);
M(1, 1) = M(1, 1) + 1 / h1 + h1 / 3;
M(n, n) = M(n, n) + 1 / hEnd + hEnd / 3;
disp([M(1,1),M(n,n)]);
end

function b = assembleLoad1D(xx,f)
n = length(xx) - 2;
b = zeros(n, 1);
fx = f(xx);
for k = 2:(n+1)
  hK = xx(k) - xx(k-1);
  hKp1 = xx(k+1) - xx(k);
  b(k-1) = fx(k) * (hK + hKp1) / 2;
end
end