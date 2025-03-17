%Main program to compute the CI to compare evenly distributed points and
%optimized points under IID noise. 
 clear;clc;
tic
 
t_best_IIDFIM={[6.34994863705261,16.7861782858528,80.4979392777307]
[6.20977450676229,15.2750903927458,18.0692952336378,80.4922051169937]
[6.27577407941707,15.7021317824451,17.7021364588466,78.4963945183739,80.4986631701891]
[5.34126363308649,7.34126604902159,15.8463055940299,17.8463085938521,78.4948727676437,80.4962916056809]
[5.33730820291159,7.33732568473259,15.8329599735807,17.8329799902337,76.4053885595053,78.4511924475796,80.4786168679362]
[5.19742621425561,7.19742714665278,14.7063616024953,16.7063635963056,18.7063650782529,76.4943377507751,78.4971540817235,80.4986760854827]
[4.34607220991320,6.34608563188164,8.34610109975754,14.9428091241092,16.9428267679213,18.9428441988758,76.4077382362805,78.4524124723076,80.4791442100305]
[4.33685461950222,6.33686782322273,8.33688291127786,14.9132214464828,16.9132380161570,18.9132539057305,74.3998089588690,76.4404162632139,78.4644728883104,80.4834737271399]
};  





% t_best_IIDSob={[8,16,80]
% [8,16,18,80]
% [8,10,16,18,80]
% [8,10,16,18,74,80]
% [8,10,16,18,74,78,80]
% [8,10,14,16,18,74,78,80]
% [8,10,12,16,18,20,74,78,80]
% [8,10,12,16,18,20,74,76,78,80]
% };

t_best_IIDSob={ [8,16,80]
[8,10,16,80]
[8,10,16,76,80]
[8,10,16,18,76,80]
[8,10,16,18,76,78,80]
[8,10,14,16,18,76,78,80]
[6,8,10,14,16,18,76,78,80]
[6,8,10,14,16,18,72,76,78,80]
};

set_numpoints=[3,4,5,6,7,8,9,10];
ave_num=1;
Sigma=0.8;%%%%consistent with OU noise Simga_OU=0.6 as phi=0.02
r=0.2; K=50; N0=4.5;
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'C_0'};   
l_set=length(t_best_IIDFIM); 
t_random=cell(l_set,1);
t_even=cell(l_set,1);
 
for i=1:l_set
    n_s=set_numpoints(i); 
    t_random{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));
    t_even{i}=linspace(t_initial, TF, n_s);%equally spaced points
end

% t_vec={t_random,t_even,t_best_IIDFIM};
t_vec={t_even,t_best_IIDFIM,t_best_IIDSob};
% prm_interest={'N_0','K'}; 
% prm_interest={'N_0'}; 
% prm_interest={'K'}; 
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
% indicating choose num_time from 30 points
% x_vec={x_OU_Sob,x_OU_FIM,x_IID_ef,x_IID_FIM};
%  t_vec={t_OU_Sob1,t_OU_Sob2,t_OU_FIM,t_IID_FIM};
% %  t_vec={t_best_OUSob,t_best_OUFIM,t_best_IIDFIM};
%   t_vec={t_best_OUSob,t_best_OUFIM,t_best_IIDFIM};


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
% Z=cell(l_set,1);
% for i=1:l_set
%       rng shuffle;
%       Z{i}=randn(ave_num,set_numpoints(i));
% end
% N_real=zeros(num_time,2);

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
%         MSE_ave_r=sqrt(mean((CI_r(2,:)-r).^2,2));
%         MSE_ave_K=sqrt(mean((CI_K(2,:)-K).^2,2));
%         MSE_ave_N0=sqrt(mean((CI_N0(2,:)-N0).^2,2));
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
%        MSE_param(k,((j-1)*num_param+1):j*num_param)=[MSE_ave_r,MSE_ave_K,MSE_ave_N0];
        suc(k,((j-1)*num_param+1):j*num_param)=[sum(r_index)/ave_num,sum(K_index)/ave_num,sum(N0_index)/ave_num];
%         MSE_est_ave(k,j)=mean(MSE_est);
    end
end
CI_parameters=CI_param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%B+FIM
%% 
% ---------------------Plot error bar 
type_name = {'Evenly','Optimized FIM','Optimized Sobol'};
num_types=length(type_name); 
lineStyles = {'*--', 's:', 'o--', 'd:'};
 
low_index = linspace(0, 0, l_set);
up_index = linspace(0, 0, l_set);
mid_index = linspace(0, 0, l_set);
CI_all=CI;
for index = 1:l_set
    low_index(index) = 4 * index - 3;
    mid_index(index) = 4 * index - 2;
    up_index(index) = 4 * index - 1;
end
middle_values=zeros(l_set,num_types);
lower_errors=middle_values;
upper_errors=middle_values;

para_dis=[0.02,4,1];
% Define color list, can be modified as per your need
colors = lines(num_types); % Using lines colormap to generate distinct colors

 fig = figure('Position', [20 20 1000 320], 'color', 'w');    
for j = 1:num_param
          subplot('Position',[0.06+(j-1)*0.31,0.13,0.265,0.75]);
%           subplot('Position', [0.06 + (j- 1) * 0.24, 0.14, 0.22, 0.75]);
             hold on;
%             ylim([CI_all(1,j)-para_dis(j),CI_all(3,j)+para_dis(j)]);    
           xlim([set_numpoints(1),set_numpoints(end)+0.5]);        
         
        % Plot additional lines on the chart
%         plot(set_numpoints', CI_param(:, (m - 1) * num_param + j), lineStyles{m}, 'MarkerSize', 8, 'LineWidth', 1.5);
          
         for i=1:num_types
             middle_values(:,i)=CI_all(mid_index, (i-1)*num_param + j);  % Middle values are in the 2nd column
             lower_errors(:,i) = middle_values(:,i) - CI_all(low_index,(i-1)* num_param + j); % Difference between middle and lower bound (column 3)
             upper_errors(:,i)= CI_all(up_index, (i-1)* num_param + j) - middle_values(:,i); % Difference between upper bound (column 1) and middle
             % Plot the error bar with enhanced style
          errorbar(set_numpoints + 0.2*(i-1), middle_values(:,i),  lower_errors(:,i),  upper_errors(:,i), lineStyles{i},...
              'MarkerFaceColor', 'c',... % Fill marker with orange
                'LineWidth', 1.5, 'MarkerSize', 8); 
%              'Color', [0.2, 0.6, 0.8],
          end
           if j == 1
%                ylabel(sprintf('Parameter %s', prm_name{j}));
               ylabel('Parameter value');
                      
               legend('Evenly','Optimized FIM','Optimized Sob','Location', 'best'); % Add legend
           end
           plot([set_numpoints,set_numpoints(end)+0.2], linspace(prm(j), prm(j), length(set_numpoints)+1),...
               'k--', 'MarkerSize', 8, 'LineWidth', 2);
        xlabel('Number of points');
        % Ensure all x-axis points are displayed
        xticks(set_numpoints);
          axis('square');  
          title(sprintf('Parameter %s', prm_name{j}));   
        grid on; hold off;
    end
   
elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
