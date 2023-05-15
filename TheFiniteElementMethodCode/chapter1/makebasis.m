function makebasis(xp,dt)
n = length(xp) - 2;
xx = xp(1):dt:xp(n+2);
rho = zeros(n+2,length(xx));
for i = 1:n
    hp = xp(i+1)-xp(i);
    hf = xp(i+2)-xp(i+1);
    
end
end