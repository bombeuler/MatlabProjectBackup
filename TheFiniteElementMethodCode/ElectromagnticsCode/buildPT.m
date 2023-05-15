function [p,t] = buildPT(insideHalfWidth,outsideHalfWidth,h)
N1 = (outsideHalfWidth - insideHalfWidth)/h;
N2 = outsideHalfWidth*2/h;
N3 = 2*insideHalfWidth/h;
Np1 = (N1+1)*(N2+1);
Np2 = 2*(N3-1)*(N1+1);
Np3 = (N1+1)*(N2+1);
p1 = zeros(2,Np1);
p2 = zeros(2,Np2);
p3 = zeros(2,Np3);
for ii = 0:N1
    for jj = 0:N2
        xx = jj*h - outsideHalfWidth;
        yy1 = ii*h - outsideHalfWidth;
        yy2 = ii*h + insideHalfWidth;
        pIndex = jj + 1 + ii*(N2+1);
        p1(:,pIndex) = [xx;yy1];
        p3(:,pIndex) = [xx;yy2];
    end
end

for ii = 1:(N3-1)
    for jj = 0:(2*N1+1)
        if(jj <= N1)
            xx = jj*h - outsideHalfWidth;
        else
            xx = (jj-N1-1)*h + insideHalfWidth;
        end
        yy = ii*h - insideHalfWidth;
        pIndex = jj + 1 + 2*(ii-1)*(N1+1);
        p2(:,pIndex) = [xx;yy];
    end
end

t1 = zeros(3,2*N1*N2);
t2 = zeros(3,4*N1*N3);
t3 = zeros(3,2*N1*N2);

for ii = 1:N1
    for jj = 1:N2
        eIndex = jj + (ii-1)*N2;
        pIndex = jj + (ii-1)*(N2+1);
        t1(:,eIndex*2 - 1) = [pIndex;pIndex+1;pIndex+N2+2];
        t1(:,eIndex*2) = [pIndex;pIndex+N2+1;pIndex+N2+2];
        p2Index = pIndex+Np1+Np2;
        t3(:,eIndex*2 - 1) = [p2Index;p2Index+1;p2Index+N2+2];
        t3(:,eIndex*2) = [p2Index;p2Index+N2+1;p2Index+N2+2];
    end
end

for ii = 1:N3
    for jj = 1:N1
        eIndex = jj +(ii-1)*N1;
        if (ii ==1)
            pIndex = jj + N1*(N2+1);
            p2Index =pIndex +N1+N3;
            t2(:,eIndex*2 - 1)=[pIndex;pIndex+1;pIndex+N2+2];
            t2(:,eIndex*2)=[pIndex;pIndex+N2+1;pIndex+N2+2];
            t2(:,2*N1*N3+eIndex*2 - 1)=[p2Index;p2Index+1;p2Index+2*N1+3];
            t2(:,2*N1*N3+eIndex*2)=[p2Index;p2Index+2*N1+2;p2Index+2*N1+3];
        elseif (ii==N3)
            pIndex = jj + 2*(ii-2)*(N1+1) + (N1+1)*(N2+1);
            p2Index =pIndex + N1+1;
            t2(:,eIndex*2 - 1)=[pIndex;pIndex+1;pIndex+2*N1+3];
            t2(:,eIndex*2)=[pIndex;pIndex+2*N1+2;pIndex+2*N1+3];
            t2(:,2*N1*N3+eIndex*2 - 1)=[p2Index;p2Index+1;p2Index+N2+2];
            t2(:,2*N1*N3+eIndex*2)=[p2Index;p2Index+N2+1;p2Index+N2+2];
        else
            pIndex = jj + 2*(ii-2)*(N1+1) + (N1+1)*(N2+1);
            p2Index =pIndex + N1+1;
            t2(:,eIndex*2 - 1)=[pIndex;pIndex+1;pIndex+2*N1+3];
            t2(:,eIndex*2)=[pIndex;pIndex+2*N1+2;pIndex+2*N1+3];
            t2(:,2*N1*N3+eIndex*2 - 1)=[p2Index;p2Index+1;p2Index+2*N1+3];
            t2(:,2*N1*N3+eIndex*2)=[p2Index;p2Index+2*N1+2;p2Index+2*N1+3];
        end
    end
end

p = [p1 p2 p3];
t = [t1 t2 t3];
end