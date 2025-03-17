reslut_x = Dec(end,:) 
Delta_t2=reslut_x(nmax+1:end-1);
N0=reslut_x(end-1);
tf=reslut_x(end)
variables=reslut_x(1:num_time)';

t1=linspace(0,tf,num_time)';
t2=t1.*variables;
 

% 使用逻辑索引来去掉0元素
t=t2(t2 ~= 0);
N_selection=fN(r,K,t,N0);
N=fN(r,K,t1,N0);
figure;
plot(t1, N, '-',t,N_selection,'o','LineWidth', 2);  
legend('model curve', 'Selection point');
points_num=sum(variables)
num_points