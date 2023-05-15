function [tt,waiting_time,inner_time,queue_length,people_numbers] = mmc(lambda,mu,n_max,c)

T_interval=exprnd(1/lambda,1,n_max);%每个人到达的间隔时间
Tc=exprnd(1/mu,1,n_max);%每个人的服务时间
c_list=zeros(1,c);%服务台状态，0代表没人在检测，数据代表对应序号的人的服务时间
queue=java.util.ArrayDeque();%初始化排队列，记录队伍各个位置的人,add()进入队尾，remove()离开队首开始检测
leave_time=zeros(1,n_max);%每个人离开队伍的时间
arrived_time=cumsum(T_interval);%每个人的到达时间

dt=0.01;
n=0;
arrived_people=0;

queue_length=zeros(1,dt*n_max);%队伍的长度
people_numbers=zeros(1,dt*n_max);%系统内的人数

while(true)
    t=n*dt;

    % 如果来人就进入队尾
    if arrived_people < n_max
        if t>=arrived_time(arrived_people+1) 
            arrived_people=arrived_people+1;
            queue.add(arrived_people);
        end
    end

    % 判断是否有空的工作台
    c_unuse=length(c_list(c_list==0));

    %有空的工作台且队列有人就进入人
    if c_unuse>0 && ~queue.isEmpty()
        place=find(c_list==0);
        l_min=min(length(place),queue.size());
        for i=1:l_min
            leave_people=queue.remove();
            c_list(place(i))=Tc(leave_people);
            leave_time(leave_people)=t;
        end
    end
    queue_length(n+1)=queue.size();
    people_numbers(n+1)=queue.size()+length(c_list(c_list~=0));

    %中止判断
    if arrived_people==n_max && queue.size()==0
        break;
    end

    %经过dt时间后服务台的变化
    c_list=max(c_list-dt,0);

    n=n+1;
end
    
waiting_time=leave_time-arrived_time;
inner_time=waiting_time+Tc;
tt=(0:n).*dt;

end