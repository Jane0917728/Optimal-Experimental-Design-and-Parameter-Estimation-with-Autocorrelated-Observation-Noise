
% Generate the data optimized using Sobol's index for Fig. 12 and Fig. 13 
r=0.2; K=50; N0=4.5;
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;% row vector

num_time=length(t);
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'N_0'};  
prm_interest={'r', 'K', 'N_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));

prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);    

percent=[0.1, 0.1, 0.1];
lower_bounds = prm'.* (1-percent);
upper_bounds = prm'.* (1+percent);

 
num_MC=20000;% number of Monte Carlo
 
S_sob=zeros(num_param,num_time);
St_sob=S_sob;
tstart=tic;
for i=1:num_time
%     SobrK(i,:)=Get_Sobol(param_list, lb_param,ub_param,N0,t(i),num_MC);
       [Total,S_sob(:,i),St_sob(:,i)]=Get_Sobol3(prm_name,lower_bounds,upper_bounds,t(i),num_MC);
end
% Sti=St_sob./max(St_sob,[],2);
% Si=S_sob;
Sti=St_sob;
% Sti=St_sob./max(St_sob,[],2);
N_value1=Nfunction(r,K,t,N0);
N_value=N_value1/max(N_value1,[],2);
figure;
hold on;
plot(t,Sti(1,:),'-o',t,Sti(2,:),'-p',t,Sti(3,:),'-d','Linewidth',2);
plot(t,N_value,':','Linewidth',2);
xlabel('time');
 xlim([0,80.5]);
ylabel('Total effect sensitivity');
% ylim([0,1]);
title ('Sobol index');
legend('r','K','C_0','C(t)');
hold off; 

figure;
hold on;

% plot(t,St_sob(1,:)/max(St_sob(2,:)),'--o',t,St_sob(2,:)/max(St_sob(2,:)),'-.p',t,St_sob(3,:)/max(St_sob(2,:)),'-d','Linewidth',1.5);
% % plot(t,N_value1/max(St_sob(2,:)),':','Linewidth',2);
% xlabel('time');
%  xlim([0,80]);
% ylabel('Total effect sensitivity');
%  ylim([0,1]);
% title ('Sobol index');
% 
% legend('r','K','C_0');
% hold off; 

 

     phi=[0.0100
    0.0206
    0.0424
    0.0874
    0.1800
    0.2300
    0.3000
    0.5250
    0.8200
    1.1150
    1.4100
    1.7050
    2.0000];

% SIG=eye(num_time,num_time)*Sigma^2;
num_run=50;
%set_numpoints=[4];
 % set_numpoints=[3,4,5,6,7,8,9,10];
numpoints=11;
num_phi=length(phi);
lbT=2; 

 
 
% x0=cell(l_set,1);


n_s=numpoints;
% lb=[t_initial,zeros(1,n_s-2)+t_initial,lbT];  % decision varaible's lower boundary
% ub=[ones(1,n_s-1)*TF,TF];


results_best=zeros(num_phi,1);
x_best=zeros(num_phi,num_time); 
count=zeros(num_phi,1);


%%%%%%%%%%%Covariance matrix for OU process 
SIG_OU=zeros(num_time,num_time); % The covariance matrix  of OU num_time*num_time 

for i=1:num_phi
    
    for s = 1:num_time
        for tt = 1:num_time
                SIG_OU(s, tt) = exp(-phi(i)*abs(t(tt)-t(s)))/(2*phi(i));  
        end   
    end 

   SIG=SIG_OU;
 
    [results_best(i),x_best(i,:)]=Opt_discret(n_s,Sti,num_time,SIG);  
    for k=1:num_run
       [results,x]=Opt_discret(n_s,Sti,num_time,SIG); 
       if results<results_best(i)
           results_best(i)=results;
           x_best(i,:)=x;
           count(i,1)=count(i,1)+1;
       end
    end
end
t_bestSob=cell(num_phi,1);
for i=1:num_phi
t_bestSob{i}=t(x_best(i,:)==1);
end