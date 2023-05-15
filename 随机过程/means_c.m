% n_max=5000;
% lambda=150;
% mu=6;
% iter_times=50;
% c=12;
% [waiting_time_mean,inner_time_mean,queue_length_mean,people_numbers_mean]=mmc_mean(n_max,lambda,mu,c,iter_times)
% disp([waiting_time_mean,inner_time_mean,queue_length_mean,people_numbers_mean]);
for ii=1:4
    saveas(ii,strcat('cdiff_',num2str(ii)),'epsc');
end
