%%%%%%%%%%%%%%%%%%Optimal experimental design for IID noise with Sob objective fucntion.
%Decision variables:  measurement points 
%This program generate the data t_best for the left-hand plots in Figure 3
%for IID noise with Sob
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

percent=[0.15, 0.15, 0.15];
lower_bounds = prm'.* (1-percent);
upper_bounds = prm'.* (1+percent); 
num_MC=20000;% number of Monte Carlo 
S_sob=zeros(num_param,num_time);
St_sob=S_sob;
% tstart=tic;
for i=1:num_time
%     SobrK(i,:)=Get_Sobol(param_list, lb_param,ub_param,N0,t(i),num_MC);
       [Total,S_sob(:,i),St_sob(:,i)]=Get_Sobol3(prm_name,lower_bounds,upper_bounds,t(i),num_MC);
end
Si=S_sob;
Sti=St_sob./max(St_sob,[],2); 
% SIG=eye(num_time,num_time)*Sigma^2;
num_run=50;
set_numpoints=[3,4,5,6,7,8,9,10];
l_set=length(set_numpoints);
t_opt=cell(l_set,1);
results_best=zeros(l_set,1);
x_best=zeros(l_set,num_time); 
count=zeros(l_set,1);

 
SIG=eye(num_time,num_time);
SIG_pinv=pinv(SIG);

for i=1:l_set
    n_s=set_numpoints(i);
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
t_bestSob=cell(l_set,1);
for i=1:l_set
t_bestSob{i}=t(x_best(i,:)==1);
end
