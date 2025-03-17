%%%%%%%%%%%%%%%%%%Optimal experimental design for IID noise with FIM objective fucntion.
%Decision variables:  measurement points 
%This program generate the data t_best for the left-hand plots in Figure 3 for IID noise with FIM
clear
clc;   %   load('norm_random.mat'); %   load('norm_rand2.mat');
r=0.2; K=50; N0= 4.5;  
% Sigma=0.8;  The variance sigma^2 does not affect the process of optimising the experimental design as it serves
% merely as a scaling factor in the objective function.
phi=0.02;

t_initial=0;
TF=80+t_initial;
min_int=2;
max_int= TF;lbT=1;
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0); 

num_run=50;
% num_soft_mod=2;
%  set_numpoints=[3,4,5,6,7,8,9,10];%set of number of points which are selected as measurement 
set_numpoints=[11];
l_set=length(set_numpoints);
%% 
count_fmin=zeros(l_set,1);
results_best=zeros(l_set,1);
result_try3=zeros(l_set,1);
t_best=cell(l_set,1);
x0=cell(l_set,1);
 
for i=1:l_set
    n_s=set_numpoints(i);
    lb=[t_initial,zeros(1,n_s-2)+t_initial,lbT];  % decision varaible's lower boundary
    ub=[ones(1,n_s-1)*TF,TF]; % decision varaible's upper boundary
    x0{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));% initial data
    [results_best(i),t_best{i}]=Opt_con_fmin_OU(r,K,N0,phi,min_int,max_int,x0{i},lb,ub);

    for k=1:num_run 
           [results,x]=Opt_con_fmin_OU(r,K,N0,phi,min_int,max_int,x0{i},lb,ub);
            x0{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));
           if results<results_best(i) 
               results_best(i)=results;               
               t_best{i}=x; 
              count_fmin(i)= count_fmin(i)+1; 
          end
    end
end

 



