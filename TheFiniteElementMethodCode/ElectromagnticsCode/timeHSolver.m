gi = [1 0 0 1 0 0 0 0 0 0];
dl = defineRegion(gi,-3,3,3,-3,2,"(PMLOUT+SQAURE)-CIRCLE",["CIRCLE","SQAURE","PMLOUT"]);
% dl = decsg(gd,"(PMLO+PMLIN)-SCAT",["PMLO","PMLIN","SCAT"]);
% dl = decsg([2 4 -5 5 5 -5 -5 -5 5 5]');
% % [p,e,t] = poimesh(dl,2,2);
[p,e,t] = initmesh(dl,"hmax",0.1);
% [p,t] = buildPT(1,5,1/6);
% t=[t;ones(1,size(t,2))]
% pdeplot(p,e,t);
% hold on
[ee,te,t2e,p2t] = buildEdges(p,t(1:3,:));
% disp(ee);
% disp(te);
np = size(p,2);
ne = size(ee,2);
nt = size(te,2);

OMEGA = 2*pi;
EPSILON = 1;
MU = 1;

A = sparse(ne,ne);
b = sparse(ne,1);

pml = buildPML(-3,3,-3,3,2,80,4);


fixed = [];
fixedin = [];
fixedinE =[];

for ek = 1:ne
    ekp = p(:,ee(:,ek));
    meanxy = mean(ekp,2);
    len = norm(ekp(:,2) - ekp(:,1));
    x = meanxy(1);
    y = meanxy(2);
    normal = -[x;y]./sqrt(x^2+y^2);
    if (x <-4.999 || x >4.999 || y <-4.999 || y >4.999)
        fixed = [fixed ek];
    elseif (sqrt(x^2+y^2) <1.001)
        fixedin = [fixedin ek];
        fixedinE = [fixedinE exp(-OMEGA*y*1i)*len];
    end


end

for kk = 1:nt
    x = p(1,te(1:3,kk));
    y = p(2,te(1:3,kk));

    SIGMA = pml(mean(x),mean(y));
% SIGMA = 0;
    KAPPA2 = 1i*OMEGA .* SIGMA - OMEGA^2 .* EPSILON; 
    [AK,bK] = buildLocalNeMat(x,y,MU,KAPPA2,[0;0]);
    A(te(4:6,kk),te(4:6,kk)) = A(te(4:6,kk),te(4:6,kk)) + AK;
    b(te(4:6,kk)) = b(te(4:6,kk)) + bK;

end

free = setdiff(1:ne,[fixed,fixedin]);
fixedAll = [fixed,fixedin];

llx = zeros(ne,1);
llx(fixedin) = fixedinE;
% llx(fixed) = exp(OMEGA*1i);
br = b(free) - A(free,fixedin) * llx(fixedin);
llx(free) = A(free,free)\br;

xxE = zeros(3*nt,1);
yyE = zeros(3*nt,1);
Exx = zeros(3*nt,1);
Eyy = zeros(3*nt,1);

xMi = zeros(nt,1);
yMi = zeros(nt,1);
kmVectorU = zeros(nt,1);
kmVectorV = zeros(nt,1);

for kk = 1:nt
%     meanxy = mean(p(:,ee(:,kk)),2);
%     xmid = meanxy(1);
%     ymid = meanxy(2);
%     areaIndex = max(t2e(:,kk));
    x = p(1,te(1:3,kk));
    y = p(2,te(1:3,kk));
    area = polyarea(x,y);
    
    xmids = mean(x);
    ymids = mean(y);
    xMi(kk) = xmids;
    yMi(kk) = ymids;
%     lens = vecnorm([x(2)-x(3) x(1)-x(3) x(1)-x(2);y(2)-y(3),x(1)-x(3) x(1)-x(2)],2);

    s1 = [y(1) - ymids;xmids - x(1)]  ./ (2*area);
    s2 = [ymids - y(2);x(2) - xmids]  ./ (2*area);
    s3 = [y(3) - ymids;xmids - x(3)]  ./ (2*area);
    smid = [s1 s2 s3] * llx(te(4:6,kk));
    kmVectorU(kk) = smid(1);
%     + exp(-1i * OMEGA * xmid);
    kmVectorV(kk) = smid(2);

    
    xmms = [x(1)/2+x(2)/4+x(3)/4 ,x(1)/4+x(2)/2+x(3)/4,x(1)/4+x(2)/4+x(3)/2];
    ymms = [y(1)/2+y(2)/4+y(3)/4 ,y(1)/4+y(2)/2+y(3)/4,y(1)/4+y(2)/4+y(3)/2];

    for mm = 1:3
        xmm = xmms(mm);
        ymm = ymms(mm);
        s1 = [y(1) - ymm;xmm - x(1)]  ./ (2*area);
        s2 = [ymm - y(2);x(2) - xmm]  ./ (2*area);
        s3 = [y(3) - ymm;xmm - x(3)]  ./ (2*area);
        smm = [s1 s2 s3] * llx(te(4:6,kk));
        Eindexs = (kk-1)*3 +mm;
        xxE(Eindexs) = xmm;
        yyE(Eindexs) = ymm;
        Exx(Eindexs) = smm(1);
        Eyy(Eindexs) = smm(2);

    end


end


figure(1);
pdegplot(dl);
hold on
pdeplot(p,[],[sort(t(1:3,:));t(4,:)],XYdata=sqrt(real(kmVectorU).^2+real(kmVectorV).^2),ColorMap="jet",XYStyle="interp",Contour="off");
hold on
% quiver(xxE,yyE,real(Exx),real(Eyy));
quiver(xMi(1:24:end),yMi(1:24:end),real(kmVectorU(1:24:end)),real(kmVectorV(1:24:end)),'w',AutoScale="on",LineWidth=0.75);
axis([-3 3 -3 3]);


figure(2);
pdegplot(dl);
hold on
pdeplot(p,[],t);
title(["nt=",num2str(nt)]);

% xMids = zeros(ne,1);
% yMids = zeros(ne,1);
% ExMids = zeros(ne,1);
% EyMids = zeros(ne,1);
% 
% for kk = 1:ne
%     meanxy = mean(p(:,ee(:,kk)),2);
%     xmid = meanxy(1);
%     ymid = meanxy(2);
%     xMids(kk) = xmid;
%     yMids(kk) = ymid;
%     Emid = zeros(2,1);
%     areaIndex = t2e(:,kk);
%     for tri = areaIndex'
%         if (tri ~=-1)
%             xp = p(1,te(1:3,tri));
%             yp = p(2,te(1:3,tri));
%             area = polyarea(xp,yp);
%             S1 = [yp(1) - ymid;xmid - xp(1)] ./ (2*area);
%             S2 = [ymid - yp(2);xp(2) - xmid] ./ (2*area);
%             S3 = [yp(3) - ymid;xmid - xp(3)] ./ (2*area);
%             Emid = Emid + [S1 S2 S3] * llx(te(4:6,tri));
%         end
%     end
%    
%     
%     ExMids(kk) = Emid(1);
%     + exp(-1i * OMEGA * xmid);
%     EyMids(kk) = Emid(2);
% end

% xxs = p(1,:);
% yys = p(2,:);
% Exs = zeros(np,1);
% Eys = zeros(np,1);
% 
% for kk = 1:np
%     xp = p(1,kk);
%     yp = p(2,kk);
% 
%     EE = zeros(2,1);
%     areaIndex = p2t(:,kk);
%     areaIndex(areaIndex == 0) = [];
%     edges = reshape(te(4:6,areaIndex),1,[]);
%     disp(areaIndex');
%     repeatEdges = setdiff(edges,unique(edges));
%     for tri = areaIndex'
%             disp(tri);
%             xx = p(1,te(1:3,tri));
%             yy = p(2,te(1:3,tri));
%             area = polyarea(xx,yy);
%             Spre = ones(3,1);
%             for ii = 1:3
%                 edge = te(3+ii,tri);
%                 if ismember(edge,repeatEdges)
%                     Spre(ii) = 1/2;
%                 end
%             end
%             S1 = Spre(1).*[yy(1) - yp;xp - xx(1)] ./ (2*area);
%             S2 = Spre(2).*[yp - yy(2);xx(2) - xp] ./ (2*area);
%             S3 = Spre(3).*[yy(3) - yp;xp - xx(3)] ./ (2*area);
%             EE = EE + [S1 S2 S3] * llx(te(4:6,tri));
%     end
%    
%     
%     Exs(kk) = EE(1);
% %     + exp(-1i * OMEGA * xmid);
%     Eys(kk) = EE(2);
% end
% 
% hh = 0.05
% 
% [Xq,Yq] = meshgrid([-5:hh:-1 1:hh:5],[-5:hh:-1 1:hh:5]);

% pdeplot(p,e,t,XYdata=real(kmVectorV),XYStyle="off",Contour="on",Levels=10);
% pdegplot(dl);
% hold on

% Vq = griddata(xxs,yys,Exs,Xq,Yq,'natural');
% 
% contour(Xq,Yq,abs(real(Vq)));



% pdesurf(p,t,ExMids);

% pdeplot(p,e,t,XYdata=abs(real(Exs)),XYStyle="off",Contour="on");
% hold on
% pdeplot(p,e,t,FlowData=[real(Exs),real(Eys)]);

% quiver(xMids,yMids,real(ExMids),real(EyMids));
% axis([-3,3,-3,3])



