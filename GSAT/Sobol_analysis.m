% Define the initial conditions
 
r=0.2; K=50;
lb_param=[0.17,45];
ub_param=[0.23,55];
lb_r=lb_param(1);
lb_K=lb_param(2);
ub_r=ub_param(1);
ub_K=ub_param(2);
num_time=20;
num_points=8;
t_initial=0.1;
tf=60;N0=1;
t=linspace(t_initial,tf,num_time)';
 

fN2_ori=@(r,K,t,N0)  (K*N0./((K-N0).*exp(-r.*t)+N0)).^2;% Define N^2(r, K)
fN_ori=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0);% Define N(r, K)
Integral_K=@(r,K,t,N0) N0*exp(r.*t).*(K+N0*(exp(r.*t)-1)-N0*(exp(r.*t)-1).*log(N0+(K-N0).*exp(-r.*t)));
Integral_r=@(r,K,t,N0) K./t.*log( K+(exp(r.*t)-1)*N0);

param_list=['r','K'];
% 使用integral2进行二重积分
% integral2的参数分别是函数句柄，x方向积分下限，x方向积分上限，y方向积分下限，y方向积分上限
Var1=zeros(num_time,1);
Var2=zeros(num_time,1);

V_r=Var1;
V_K=Var1;
for i=1:num_time
    fN2=@(r,K)  fN2_ori(r,K,t(i),N0);% Define N^2(r, K) for integral
    fN=@(r,K)  fN_ori(r,K,t(i),N0);% Define N(r, K) for integral
    E_K=@(r) ( Integral_K(r,ub_K,t(i),N0)-Integral_K(r,lb_K,t(i),N0)).^2./(ub_K-lb_K)^2;
    E_r=@(K) ( Integral_r(ub_r,K,t(i),N0)-Integral_r(lb_r,K,t(i),N0)).^2./(ub_r-lb_r)^2;

    Var1(i) = 1/((ub_param(2)-lb_param(2))*(ub_param(1)-lb_param(1)))*integral2(fN2, lb_param(1), ub_param(1), lb_param(2), ub_param(2));
    Var2(i) = (1/((ub_param(2)-lb_param(2))*(ub_param(1)-lb_param(1)))*integral2(fN, lb_param(1), ub_param(1), lb_param(2), ub_param(2))).^2;
    
      
        V_K(i)= 1/(ub_K-lb_K)  * integral(E_r,lb_param(2),ub_param(2)) ;
    V_r(i)= 1/(ub_r-lb_r)  * integral(E_K,lb_param(1),ub_param(1)) ;
end


V_total =Var1-Var2;
Vr=V_r-Var2;
VK=V_K-Var2;
S_r=(V_r-Var2)./V_total; %First-order 
S_K=(V_K-Var2)./V_total;
S=[S_r,S_K];
% 显示结果
% disp(['二重积分的结果是：', num2str(result)]);
