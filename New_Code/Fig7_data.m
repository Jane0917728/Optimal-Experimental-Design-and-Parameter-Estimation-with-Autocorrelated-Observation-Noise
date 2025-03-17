clear
clc;
%   load('norm_random.mat');
%   load('norm_rand2.mat');
r=0.2;
K=50;
N0= 4.5; 
t_initial=0;
TF=80+t_initial;
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0); 

num_run=5; 
numpoints=5; 
num_phi=30;
phi=exp(linspace(log(0.01),log(2.5),num_phi))';  
min_int=1;
max_int= TF;
lbT=1; 
 
% x0=cell(l_set,1);
count_fmin=zeros(num_phi,1);
results_best=zeros(num_phi,1);

t_best=cell(num_phi,1);
% Constraints=cell(l_set,num_phi);
   n_s=numpoints;
    lb=[t_initial,zeros(1,n_s-2)+t_initial,lbT];  % decision varaible's lower boundary
    ub=[ones(1,n_s-1)*TF,TF];      
    for j=1:num_phi           
%           x0{i}=linspace(t_initial,TF,n_s); 
          x0=sort(t_initial+(TF-t_initial)*rand(1,n_s)); 
         [results_best(j,1),t_best{j,1}]=Opt_con_fmin(r,K,N0,phi(j),min_int,max_int, x0,lb,ub);
 
        for k=1:num_run     
              x0=sort(t_initial+(TF-t_initial)*rand(1,n_s));
                      [results,x]=Opt_con_fmin(r,K,N0,phi(j),min_int,max_int, x0,lb,ub);
                            if results<results_best(j,1) 
                                results_best(j,1) =results;
                                t_best{j,1}=x;
                                count_fmin(j,1)= count_fmin(j,1)+1;
                                x0=x; 
                            else
                                 x0=sort(t_initial+(TF-t_initial)*rand(1,n_s));
                            end 
        end
    end 

 