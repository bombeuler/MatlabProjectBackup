n_max=5000;
lambda=150;
mu=6;
c=12;

[tt,waiting_time,inner_time,queue_length,people_numbers] =mmc(lambda,mu,n_max,c);

figure(1);
plot(tt,queue_length,'-b');
title('队列长度');
xlabel('时间（单位：min）');
ylabel('人数');
figure(2);
plot(tt,people_numbers,'-r');
title('系统内人数');
xlabel('时间（单位：min）');
ylabel('人数');
figure(3);
plot(waiting_time,"-g");
title('等待时间');
xlabel('人的序号');
ylabel('时间（单位：min）');
figure(4);
plot(inner_time,"-k");
title('停留时间');
xlabel('人的序号');
ylabel('时间（单位：min）');

for ii=1:4
    saveas(ii,strcat('plot_mmc_',num2str(ii)),'epsc');
end