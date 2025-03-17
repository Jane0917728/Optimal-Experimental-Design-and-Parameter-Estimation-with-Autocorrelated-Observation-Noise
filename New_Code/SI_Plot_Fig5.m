 clear;clc;
tic
set_numpoints=[5,6,7,8,9];
ave_num=1;
Sigma=0.8;%%%%consistent with OU noise Simga_OU=0.6 as phi=0.02
r=0.2; K=50; N0=4.5;
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'C_0'};   
l_set=length(set_numpoints); 
t_random=cell(l_set,1);
t_even=cell(l_set,1);
 
for i=1:l_set
    n_s=set_numpoints(i); 
%     t_random{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));
    t_even{i}=linspace(t_initial, TF, n_s);%equally spaced points
end
% t_check={[0,20,40,60,80]
% [0,20,40,48,64,80]
% [0,13.3333333333333,20,40,66.6666666666667,80]
% [0,11.4285714285714,20.8571428571429,40.2857142857143,45.7142857142857,57.1428571428571,68.5714285714286,80]
% [0,10,20,30,40,50,60,70,80]
% };

% t_check={[0,20,40,60,80]
% [0,20,40,48,64,80]
% [0,13.3333333333333,20,40,66.6666666666667,80]
% [0,11.4285714285714,20.8571428571429,40.2857142857143,45.7142857142857,57.1428571428571,68.5714285714286,80]
% [0,10,20,30,40,50,60,70,80]
% };

t_check={[0,20,40,60,80]
[0,20,40,60,64,80]
[0,13.3,20,40,60,80]
[0,11.4,20,40,45.7,60,68.6,80]
[0,10,20,30,40,50,60,70,80]
};
% t_vec={t_random,t_even,t_best_IIDFIM};
t_vec={t_even,t_check};

prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

%%%%%%%%%%%parameters for OU process noise
lb_p=[0.05,25,0.1]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.4,75,15];

CI_bounds=[0.015,1,1;
                        0.065,3.7,3.5];

numpts=70; % parameter number of the profile likelihood
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  
symbols = ['o', 'p', '+', '^', 'x', 's', 'd',  'v', '<', '>', 'h'];



% x_vec={x_vecBS,x_vecBF,x_vecIIDSob,x_vecIIDFIM};
num_types=length(t_vec); %number of types, such as OU noise with FIM, or OU noise with Sob index.
CI_w=zeros(l_set,num_param*num_types);%The first  means the number of change parameters
CI=zeros(l_set*4,num_param*num_types); %4 denotes all information CI_lower bound, estimated value, 
%CI upper bound and width of CI.
 
CI_width=zeros(l_set,num_param*num_types);
CI_param=CI_width;
MSE_param=CI_w; % MSE of each parameters
suc=zeros(l_set,num_param*num_types); %8*4
MSE_est_ave=zeros(l_set,num_param*num_types); %the MSE of the estimated parameters
CI_r=zeros(3,ave_num);  %CI_r denotes the lower bound, the estimated value and the upper bound.
CI_K=CI_r; % CI_K denotes the lower bound, the estimated value and the upper bound.
CI_N0=CI_r; % CI_N0 denotes the lower bound, the estimated value and the upper bound.
 CI_1=CI;
 CI_w1=CI_w;
MSE_est=zeros(ave_num,1); % record the mean of square error between the estimated value and the real value.


for j=1:num_types
    t_v=t_vec{j};
    for k=1:l_set
        t_opt=t_v{k,1};
        N_real=Nfunction(r,K,t_opt,N0);
        n_s=length(t_opt);
        
         SIG=eye(n_s,n_s)*Sigma^2; % IID covariance matrix             
 
        parfor i=1:ave_num    
               measurements=N_real+ Sigma*randn(size(N_real)); %  IID measurement noise
             [CI_r(:,i),CI_K(:,i),CI_N0(:,i),MSE_est(i)] = ...
                 get_CI(r,K,N0,numpts,SIG,lb_p,ub_p,measurements,t_opt,num_param);
        end            
  
         r_index=(CI_r(1,:)>-0.2)&(CI_r(3,:)<1);
        K_index=(CI_K(1,:)>10)&(CI_K(3,:)<110); 
        N0_index=(CI_N0(1,:)>-1)&(CI_N0(3,:)<20);
         %%%%%%%%%%%%%%%%%%%%%%%% average for all the results of ave_num
        %%%%%%%%%%%%%%%%%%%%%%%% running
        CI_r_ave1 = mean(CI_r,2);
        CI_K_ave1= mean(CI_K,2);
        CI_N0_ave1 = mean(CI_N0,2);
        
        CI_r_ave = mean(CI_r(:,r_index),2);
        CI_K_ave= mean(CI_K(:,K_index),2);
        CI_N0_ave = mean(CI_N0(:,N0_index),2);

         %%%%%%%%%%%%%%%%%%%%%%%% average only for the results of
         % successfully getting the confident intervel.
%         CI_r_ave = mean(CI_r(:,(CI_r(1,:)>-0.2)&(CI_r(3,:)<1)),2);
%         CI_K_ave = mean(CI_K(:,(CI_K(1,:)>10)&(CI_K(3,:)<110)),2);         
        CI_f_r1=[CI_r_ave1;CI_r_ave1(3)-CI_r_ave1(1)];
        CI_f_K1=[CI_K_ave1;CI_K_ave1(3)-CI_K_ave1(1)];
        CI_f_N01=[CI_N0_ave1;CI_N0_ave1(3)-CI_N0_ave1(1)];
           CI_f_r=[CI_r_ave;CI_r_ave(3)-CI_r_ave(1)];
        CI_f_K=[CI_K_ave;CI_K_ave(3)-CI_K_ave(1)];
        CI_f_N0=[CI_N0_ave;CI_N0_ave(3)-CI_N0_ave(1)];
        
        CI_1(4*(k-1)+1:4*k,((j-1)*num_param+1):j*num_param)=[CI_f_r1,CI_f_K1,CI_f_N01];%CI
        % CI_width and CI_param
        CI_w1(k,((j-1)*num_param+1):j*num_param)=CI_1(4*k,((j-1)*num_param+1):j*num_param);   
        CI(4*(k-1)+1:4*k,((j-1)*num_param+1):j*num_param)=[CI_f_r,CI_f_K,CI_f_N0];%CI
        % CI_width and CI_param
        CI_w(k,((j-1)*num_param+1):j*num_param)=CI(4*k,((j-1)*num_param+1):j*num_param);    
        
        CI_param(k,((j-1)*num_param+1):j*num_param)=CI(4*(k-1)+2,((j-1)*num_param+1):j*num_param);  
        
       % compute the Identification success ratio
        suc(k,((j-1)*num_param+1):j*num_param)=[sum(r_index)/ave_num,sum(K_index)/ave_num,sum(N0_index)/ave_num];

    end
end
 CI_parameters=CI_param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%B+FIM
%plot the figures
 
fig=figure('Position',[20 20 1000 320],'color','w');
  %   [left bottom width height]
for j=1:num_param
%       subplot('Position',[0.06+(j-1)*0.31,0.12,0.265,0.75]);
       subplot('Position',[0.06+(j-1)*0.31,0.13,0.265,0.75]);
%     subplot(1,num_param,j,'Position',[0.05+(j-1)*0.3,0.12,0.26,0.75]);
 %                     set_numpoints, CI_w(:,j+num_param*(num_types-1)), 'o--', ...
        plot(set_numpoints(1:end),  CI_w(1:end,j+0),'*-',set_numpoints,CI_w(:,j+num_param*(num_types-1)), 'x:',...
                                                 'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);

    
    xlim([set_numpoints(1),set_numpoints(end)]);
     ylim([CI_bounds(1,j),CI_bounds(2,j)]);
    xlabel('Number of Points');
    ylabel('Confifence Intervals Width');
%        legend('IID FIM','OU FIM','OU Sobol','Location','Northeast');
%     legend('Evenly','Test1','Test2'); 
       legend('Evenly','Test'); 
    title(sprintf('Parameter %s', prm_name{j}));   
end

   
elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
