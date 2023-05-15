dl=Scatterg();
[p,dummy,t]=initmesh(dl,'hmax',0.1);
[ee,te,t2e,p2t] = buildEdges(p,t(1:3,:));
neighbors = Tri2Tri(p,t);
e = Tri2Edge(p,t);
ne = max(e(:));
nt = size(t,2);
np = size(p,2);
A = sparse(ne,ne);
b = zeros(ne,1);
omega = 2*pi/1;
mu = 1;
epsilon = 1;
signs=2*((neighbors'<[1:nt; 1:nt; 1:nt])-1/2)
for i=1:size(t,2)
    nodes = t(1:3,i);
    x = p(1,nodes);
    y = p(2,nodes);
    edges = e(i,:);
    sigma = 0;
    sd=t(4,i);
    if sd==7 % cavity
        %
    elseif sd==2 | sd==9 % up down pml
        sigma=5*(abs(mean(y))-3)^4;
    elseif sd==6 | sd==8 % left right pml
        sigma=5*(abs(mean(x))-3)^4;
    else
        sigma=5*(abs(mean(x))-3)^4+(abs(mean(y))-3)^4;
    end
    kappa=sqrt(sqrt(-1)*sigma*omega-epsilon*omega^2);
    [AE,FE] = RotRot(x,y,signs(:,i),mu,kappa,[0,0]);
    A(edges,edges) = A(edges,edges) + AE;
    b(edges) = b(edges) + FE;
end
[xmid,ymid,edges] = EdgeMidPoints(p,e,t);
fixed=[]; % fixed nodes
gvals=[];
for i=1:ne % loop over edges
    r=edges(i); % edge or node number
    x=xmid(i); % node x-coordinate
    y=ymid(i); %
    if (sqrt(x^2+y^2)<1.001) % cylinder
        normal = -[x; y]/sqrt(x^2+y^2);
        fixed = [fixed; r];
        gvals = [gvals; normal(2)*exp(omega*y*sqrt(-1))];
    end
    if (abs(x)>4.999 | abs(y)>4.999) % pml
        fixed = [fixed; r];
        gvals = [gvals; 0];
    end
end
free=setdiff([1:ne],fixed);
b=b(free)-A(free,fixed)*gvals;
A=A(free,free);
xi=zeros(ne,1);
xi(fixed)=gvals;
xi(free)=A\b;
disp(xi);

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
%     repeatEdges = setdiff(edges,unique(edges));
%     for tri = areaIndex'
%             disp(tri);
%             xx = p(1,te(1:3,tri));
%             yy = p(2,te(1:3,tri));
%             area = polyarea(xx,yy);
%             Spre = ones(3,1);
%             for ii = 1:3
%                 edge = te(3+ii,tri);
%                 if ismember(edge,edges)
%                     Spre(ii) = 1/2;
%                 end
%             end
%             S1 = Spre(1).*[yy(1) - yp;xp - xx(1)] ./ (2*area);
%             S2 = Spre(2).*[yp - yy(2);xx(2) - xp] ./ (2*area);
%             S3 = Spre(3).*[yy(3) - yp;xp - xx(3)] ./ (2*area);
%             EE = EE + [S1 S2 S3] * xi(te(4:6,tri));
%     end
%    
%     
%     Exs(kk) = EE(1);
%     + exp(-1i * OMEGA * xmid);
%     Eys(kk) = EE(2);
% end
% 
% 
% 
% 
% pdeplot(p,e,t,XYdata=real(kmVectorV),XYStyle="off",Contour="on",Levels=10);
% pdegplot(dl);
% hold on
% 
% pdeplot(p,e,t,XYdata=abs(real(Exs)),XYStyle="off",Contour="on");
% hold on
% pdeplot(p,e,t,FlowData=[real(Exs),real(Eys)]);


kmVectorU = zeros(nt,1);
kmVectorV = zeros(nt,1);

for kk = 1:nt
%     meanxy = mean(p(:,ee(:,kk)),2);
%     xmid = meanxy(1);
%     ymid = meanxy(2);
%     areaIndex = max(t2e(:,kk));
    x = p(1,t(1:3,kk));
    y = p(2,t(1:3,kk));
    area = polyarea(x,y);
    edges = e(kk,:);
    signss = signs(:,kk);
    xmids = mean(x);
    ymids = mean(y);
    lens = vecnorm([x(2)-x(3) x(1)-x(3) x(1)-x(2);y(2)-y(3),x(1)-x(3) x(1)-x(2)],2);

    s1 = [y(1) - ymids;xmids - x(1)] .* signss(1) .* lens(1) ./ (2*area);
    s2 = [ymids - y(2);x(2) - xmids] .* signss(2) .* lens(2) ./ (2*area);
    s3 = [y(3) - ymids;xmids - x(3)] .* signss(3) .* lens(3) ./ (2*area);
    smid = [s1 s2 s3] * xi(e(kk,:));
    kmVectorU(kk) = smid(1);
%     + exp(-1i * OMEGA * xmid);
    kmVectorV(kk) = smid(2);
end

pdegplot(dl);
hold on
pdeplot(p,e,t,XYdata=real(kmVectorU),XYStyle="off",Contour="on",Levels=10);
hold on
% pdeplot(p,e,t,XYdata=real(xi),XYStyle="off",Contour="on");
% 
% % [xq,yq] = meshgrid(-5:.5:5, -5:.5:5);
% % vq = griddata(xmid,ymid,real(xi),xq,yq);
% % 
% % contour(xq,yq,vq);
% 
pdeplot(p,e,t,FlowData=[real(kmVectorU),real(kmVectorV)]);
