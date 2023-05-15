hScale = 1.2;

figure(1);
hold on
xlabel("hmax");
ylabel("cond");
legend("on");

condMat = zeros(6,10);
hh = zeros(1,10);
for jj=0:5
    hMax = 0.2;
    for ii=1:10
        hh(ii) = hMax;
        model = createpde(2);
        geometryFromEdges(model,@lshapeg);
        meshes = generateMesh(model,"hmax",hMax,"GeometricOrder","linear","Hedge",{[1,2,3,4,5,6,9,10],hMax/(2^jj)});
        p = meshes.Nodes;
        t = meshes.Elements;
        M = MassAssembler2D(p, t);
        condMat(jj+1,ii) = cond(M);
    % pdeplot(meshes, ...
    %         'FaceAlpha', 0.4, 'ColorMap','jet','Mesh','on');
%     plot(hMax,cond(M),"ob");
%     hold on
        hMax = hMax / hScale;
    end
end

for jj = 1:6
    plot(hh,condMat(jj,:),"DisplayName",num2str(2^(jj-1)));
    hold on
end