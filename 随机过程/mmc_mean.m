function [waiting_time_mean,inner_time_mean,queue_length_mean,people_numbers_mean]=mmc_mean(n_max,lambda,mu,c,iter_times)
waiting_time_mean=0;
inner_time_mean=0;
queue_length_mean=0;
people_numbers_mean=0;
for simu =1:iter_times
    [tt,waiting_time,inner_time,queue_length,people_numbers] =mmc(lambda,mu,n_max,c);
    waiting_time_mean=mean([waiting_time_mean,waiting_time]);
    inner_time_mean=mean([inner_time_mean,inner_time]);
    queue_length_mean=mean([queue_length_mean,queue_length]);
    people_numbers_mean=mean([people_numbers_mean,people_numbers]);
end
end