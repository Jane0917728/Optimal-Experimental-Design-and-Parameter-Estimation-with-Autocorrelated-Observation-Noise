%test get_Sobol
clear; 
clc;
 
r=0.2; K=50;
N0=4.5;
percent=0.3;
lb_param=[r*(1-percent),K*(1-percent),N0*(1-percent)];
ub_param=[r*(1+percent),K*(1+percent),N0*(1+percent)];

num_points=8;
t_interval=2.5;
t_initial=0.1;
TF=80;
t=t_initial:t_interval:ceil(TF/t_interval)*t_interval+t_initial;
t=t';
num_time=length(t);
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0);
 
param_list={'r','K','N0'};
num_param=length(param_list);
num_MC=20000;% number of Monte Carlo
%  measurements=f(r,K,t,N0)+epsilon;
% SobrK=zeros(num_time,4);
% for i=1:num_time
%     SobrK(i,:)=Get_Sobol(param_list, lb_param,ub_param,N0,t(i),num_MC);
% end
tstart=tic;
for i=1:num_time
%     SobrK(i,:)=Get_Sobol(param_list, lb_param,ub_param,N0,t(i),num_MC);
       [S_sob(:,i),St_sob(:,i)]=Get_Sobol3(param_list,lb_param,ub_param,t(i),num_MC);
end
 elapsed1=toc(tstart);
% get the sensitivity from eFAST method
tstart=tic;
 [S_efast,St_efast]=get_efast(t);
 elapsed2=toc(tstart);
 
 %get the contribution of each time points of FIM
pNpK=@(r,K,t,N0) N0^2*(1-exp(-r*t))./((K-N0)*exp(-r*t)+N0).^2;
pNpr=@(r,K,t,N0) K*N0*(K-N0)*exp(-r*t).*t./((K-N0)*exp(-r*t)+N0).^2;
pNpN0=@(r,K,t,N0) K^2*exp(-r*t)./((K-N0)*exp(-r*t)+N0).^2;
d_param=zeros(num_time,3);
param=zeros(1,3);
d_param(:,1)=pNpr(r,K,t,N0)*r;
d_param(:,2)=pNpK(r,K,t,N0)*K;
d_param(:,3)=pNpN0(r,K,t,N0)*N0;

S_FIM=d_param';



  N_selection=fN(r,K,t,N0);
% for i=1:num_param
%   figure;
%   plot(t, St_sob(i,:)/max(St_sob(i,:)), '*',t, St_efast(i,:)/max(St_efast(i,:)),'<',t, S_FIM(i,:)/max(S_FIM(i,:)),'d','LineWidth', 2);  
% %plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<','LineWidth', 2); 
%   xlabel('time');
%   ylabel('Sensitivity');
%   legend('Sob', 'eFast', 'FIM');
% % legend('Sob', 'eFast');
% % title('IID noise');
%   title(param_list{i});
%   ax = gca; %  
%   ax.Position = [0.1 0.1 0.85 0.825]; 
% end

  figure;
 
  plot(t,  St_efast(1,:),'-<', t,  St_efast(2,:), '-d', t,  St_efast(3,:),'-o','LineWidth', 2);  
%plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<','LineWidth', 2); 
  xlabel('time');
  ylabel('Total effect sensitivity');
  legend('r', 'K', 'N0');
  title('Sensitivity using eFAST');
  
  
  figure;
 
  plot(t,  S_FIM(1,:),'-<',t,  S_FIM(2,:),'-d',t,  S_FIM(3,:),'-o','LineWidth', 2);  
%plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<','LineWidth', 2); 
  xlabel('time');
  ylabel('Total effect sensitivity');
  legend('r', 'K', 'N0');
  title('Sensitivity using FIM');
 
  
  figure;
    plot(t,  St_sob(1,:),'-<',t,  St_sob(2,:),'-d',t,  St_sob(3,:),'-o','LineWidth', 2);  
%plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<','LineWidth', 2); 
  xlabel('time');
  ylabel('Total effect sensitivity');
  legend('r', 'K', 'N0');
  title('Sensitivity using Sob');
  
  
% for i=1:num_param
% figure;
%  plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<',t, S_FIM(i,:),'d','LineWidth', 2);  
% plot(t, St_sob(i,:), '*',t, St_efast(i,:),'<','LineWidth', 2); 
% xlabel('time');
% ylabel('Sensitivity');
% legend('Sob', 'eFast', 'FIM');
% % legend('Sob', 'eFast');
% % title('IID noise');
%   title(param_list(i));
%    ax = gca; %  
% ax.Position = [0.1 0.1 0.85 0.825]; 
% end
% figure;
%   N_selection=fN(r,K,t,N0);
%  plot(t,N_selection,'-', 'LineWidth',2)
 
 