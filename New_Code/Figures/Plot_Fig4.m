clear all;
close all;

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

addpath('../GSAT/');

r=0.2;
K=50;
N0=4.5;
prm_all=[r K N0]';
prm_name_all = {'r', 'K', 'N_0'};
prm_interest={'r', 'K', 'N_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

phi=0.02;
t_interval=1;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:TF;
num_time=length(t);

% Calculate the FIM row
pNpK=@(r,K,t,N0) N0^2*(1-exp(-r*t))./((K-N0)*exp(-r*t)+N0).^2;% partial N/partial K
pNpr=@(r,K,t,N0) K*N0*(K-N0)*exp(-r*t).*t./((K-N0)*exp(-r*t)+N0).^2; % partial N/partial r
pNpN0=@(r,K,t,N0) K^2*exp(-r*t)./((K-N0)*exp(-r*t)+N0).^2; % partial N/partial N_0

d_param=zeros(num_param,num_time);
param=zeros(3,1);
param(1)=r;
param(2)=K;
param(3)=N0;
d_param(1,:)=pNpr(r,K,t,N0)*param(1);
d_param(2,:)=pNpK(r,K,t,N0)*param(2);
d_param(3,:)=pNpN0(r,K,t,N0)*param(3);
S_FIM=d_param;

figure(1);
hFig = gcf;
set(hFig,'Position',[20 20 1000 320],'color','w');
hold on;

k=1;
subplot('Position',[0.06+(k-1)*0.31,0.14,0.265,0.75]);
plot(t,S_FIM(1,:),'-',t,S_FIM(2,:),'--',t,S_FIM(3,:),'-.','Linewidth',2.0);
xlabel('Time');
ylabel('FIM sensitivity')
% legend('r','K','C_0','Location','southeast');
ylim([0 51])
xticks([0 20 40 60 80]);
yticks([0 10 20 30 40 50]);

% Plot Sobol's index for optimization
percent=[0.15, 0.15, 0.15];
lower_bounds = prm'.* (1-percent);
upper_bounds = prm'.* (1+percent);
num_MC=20000; % number of Monte Carlo
S_sob=zeros(num_param,num_time);
St_sob=S_sob;
tstart=tic;
for i=1:num_time
    [Total,S_sob(:,i),St_sob(:,i)]=Get_Sobol3(prm_name,lower_bounds,upper_bounds,t(i),num_MC);
end

k=2;
subplot('Position',[0.06+(k-1)*0.35,0.14,0.265,0.75]);
plot(t,St_sob(1,:)/max(St_sob(2,:)),'-',...
    t,St_sob(2,:)/max(St_sob(2,:)),'--',...
    t,St_sob(3,:)/max(St_sob(2,:)),'-.','Linewidth',2.0);
xlabel('Time');
ylabel('Total effect sensitivity');
ylim([0 1.02])
xticks([0 20 40 60 80]);
yticks([0 0.2 0.4 0.6 0.8 1.0]);
legend('r','K','C_0','Location','southeast');
ytickformat('%.1f');

print('Fig4','-dpng','-r300')





