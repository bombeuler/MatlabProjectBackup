function xsMat = makeMatrix()
sList = 0:15;
sValues = [0,0,5,5,10,10,15,15,20,20,25,25,30,30,35,35];
xList = 0:5;
xWeight = 3;
xValues = [0,7,14,21,28,35];
sLength = length(sList);
xLength = length(xList);
xsMat = -1 .* ones(xLength,sLength);

for sI = 1:sLength
    for xI = 1:xLength
        sIndex = sI - xWeight * xList(xI);
        if sIndex < 1
            continue;
        end

        xsMat(xI,sI) = xValues(xI) + sValues(sIndex);
    end
end

end