n_max=500;
lambda=150;
mu=6;
iter_times=10;

cc=12:20;
lcc=length(cc);
wtm_list=zeros(1,lcc);
itm_list=zeros(1,lcc);
qlm_list=zeros(1,lcc);
pnm_list=zeros(1,lcc);

for c=cc
    [wtm,itm,qlm,pnm]=mmc_mean(n_max,lambda,mu,c,iter_times);
    ic=c-11;
    wtm_list(ic)=wtm;
    itm_list(ic)=itm;
    qlm_list(ic)=qlm;
    pnm_list(ic)=pnm;
end

figure(1);
plot(cc,wtm_list,'ob');
title('平均等待时间');
xlabel('c');
ylabel('时间（单位：min）');
figure(2);
plot(cc,itm_list,'or');
title('平均停留时间');
xlabel('c');
ylabel('时间（单位：min）');
figure(3);
plot(cc,qlm_list,'og');
title('队伍平均长度');
xlabel('c');
ylabel('人数');
figure(4);
plot(cc,pnm_list,'ok');
title('系统平均人数');
xlabel('c');
ylabel('人数');
for ii=cc
    saveas(ii-11,strcat('cdiff_',num2str(ii)),'epsc');
end
