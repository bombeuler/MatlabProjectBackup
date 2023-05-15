function [ee,te,t2e,p2t] = buildEdges(p,tmod)
nt = size(tmod,2);
np = size(p,2);
eMat = sparse(np,np);
%     tii = sort(tmod(:,ii));
for ii = 1:nt
   eMat(tmod(1,ii),tmod(2,ii)) = 1;
   eMat(tmod(2,ii),tmod(3,ii)) = 1;
   eMat(tmod(3,ii),tmod(1,ii)) = 1;
end

[ei,ej,ev] = find(triu(eMat+eMat'));
ne = length(ei);
ee = [ei,ej]';
te = zeros(6,nt);
eMatNew = sparse(ei,ej,1:ne,np,np);
for ii = 1:nt
    tii = sort(tmod(:,ii));
    te(1:3,ii) = sort(tmod(:,ii));
    te(4:6,ii) = [eMatNew(tii(2),tii(3));eMatNew(tii(1),tii(3));eMatNew(tii(1),tii(2))];
end

t2e = -ones(2,ne);

for ii = 1:nt
    for jj = 4:6
        eindex = te(jj,ii);
        if(t2e(1,eindex) > 0)
            t2e(2,eindex) = ii;
        else
             t2e(1,eindex) = ii;
        end
    end
end

eMatFull = eMatNew + eMatNew';

p2t = zeros(10,np);

for ii = 1:np
    [eii,ejj,evv] = find(eMatFull(ii,:));
    tIndexs = unique(reshape(t2e(:,evv),[],1));
    if tIndexs(1) == -1
        tIndexs = tIndexs(2:end);
    end
    len = length(tIndexs);
    p2t(1:len,ii) = tIndexs;
end

end


% disp(ee);
% disp(ev);
% disp([length(ei),nt,np,size(eMat)]);

