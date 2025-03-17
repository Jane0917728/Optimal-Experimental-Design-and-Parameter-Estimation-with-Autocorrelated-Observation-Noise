% clear;
% clc; 
r=0.2; K=50; N0=4.5;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'N_0'};  
prm_interest={'r', 'K', 'N_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);    


phi=0.02;% parameter of OU
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:TF;% row vector 
num_time=length(t);


%Calculate the FIM row
pNpK=@(r,K,t,N0) N0^2*(1-exp(-r*t))./((K-N0)*exp(-r*t)+N0).^2;% partial N/partial K 
pNpr=@(r,K,t,N0) K*N0*(K-N0)*exp(-r*t).*t./((K-N0)*exp(-r*t)+N0).^2; % partial N/partial r 
pNpN0=@(r,K,t,N0) K^2*exp(-r*t)./((K-N0)*exp(-r*t)+N0).^2; % partial N/partial N_0
% pNpt=@(r,K,t,N0) N0*K*r*(K-N0)*exp(-r*t)./((K-N0)*exp(-r*t)+N0).^2;

d_param=zeros(num_param,num_time);
param=zeros(3,1);
param(1)=r;
param(2)=K;
param(3)=N0;
% dNdt=pNpt(r,K,t,N0);
d_param(1,:)=pNpr(r,K,t,N0)*param(1);
d_param(2,:)=pNpK(r,K,t,N0)*param(2);
d_param(3,:)=pNpN0(r,K,t,N0)*param(3);
 

S_FIM=d_param;
figure;
hold on;

plot(t,S_FIM(1,:),'--o',t,S_FIM(2,:),'-.p',t,S_FIM(3,:),'-d','Linewidth',1.5);

xlabel('time');
ylabel('FIM Sensitiity')
legend('r','K','C_0');
title ('Local Sensitivity based on FIM')
hold off;

%Plot Sobol's index for optimization
percent=[0.15, 0.15, 0.15];
lower_bounds = prm'.* (1-percent);
upper_bounds = prm'.* (1+percent);
num_MC=20000;% number of Monte Carlo
S_sob=zeros(num_param,num_time);
St_sob=S_sob;
tstart=tic;
for i=1:num_time
       [Total,S_sob(:,i),St_sob(:,i)]=Get_Sobol3(prm_name,lower_bounds,upper_bounds,t(i),num_MC);
end
% Sti=St_sob./max(St_sob,[],2);
% Sti=St_sob;
figure;
hold on;
%normalized;
plot(t,St_sob(1,:)/max(St_sob(2,:)),'--o',t,St_sob(2,:)/max(St_sob(2,:)),'-.p',t,St_sob(3,:)/max(St_sob(2,:)),'-d','Linewidth',1.5);
% plot(t,Sti(1,:),'-o',t,Sti(2,:),'-p',t,Sti(3,:),'-d','MarkerSize',4,'Linewidth',1.5);
xlabel('time');
ylabel('Total effect sensitivity');

title ('Sobol index');

legend('r','K','C_0');
hold off; 





