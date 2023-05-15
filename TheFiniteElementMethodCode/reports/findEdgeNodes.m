function edges = findEdgeNodes(p)
edges = [];
np = size(p,2);
for kp = 1:np
    x = p(1,kp);
    y = p(2,kp);
    if (y == 0 && x >= -1 && x <= 0)|| x == 1 || x == -1 || y == 1
        edges = [edges kp];
    end
end
end

