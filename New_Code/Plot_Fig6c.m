%Main program to compute the CI to compare evenly distributed points and
%optimized points under misspecified noise model
clear;clc;
tic
 
t_best_OUSob={[10,16,32]
[0,10,16,32]
[0,10,16,28,80]
[0,8,12,16,28,80]
[0,8,10,14,18,30,80]
[0,8,10,14,18,26,34,80]
[0,4,8,10,14,18,26,34,80]
[0,4,8,10,14,16,18,26,36,80]};%%%%%%%%%%%%

t_best_IIDSob={[8,16,80]
[8,10,16,80]
[8,10,16,78,80]
[8,10,16,18,78,80]
[8,10,16,18,74,78,80]
[8,10,14,16,18,74,78,80]
[8,10,12,16,18,20,74,78,80]
[8,10,12,16,18,20,72,74,78,80]
};
 

t_best_IIDFIM={[6.34994863705261,16.7861782858528,80.4979392777307]
[6.20977450676229,15.2750903927458,18.0692952336378,80.4922051169937]
[6.27577407941707,15.7021317824451,17.7021364588466,78.4963945183739,80.4986631701891]
[5.34126363308649,7.34126604902159,15.8463055940299,17.8463085938521,78.4948727676437,80.4962916056809]
[5.33730820291159,7.33732568473259,15.8329599735807,17.8329799902337,76.4053885595053,78.4511924475796,80.4786168679362]
[5.19742621425561,7.19742714665278,14.7063616024953,16.7063635963056,18.7063650782529,76.4943377507751,78.4971540817235,80.4986760854827]
[4.34607220991320,6.34608563188164,8.34610109975754,14.9428091241092,16.9428267679213,18.9428441988758,76.4077382362805,78.4524124723076,80.4791442100305]
[4.33685461950222,6.33686782322273,8.33688291127786,14.9132214464828,16.9132380161570,18.9132539057305,74.3998089588690,76.4404162632139,78.4644728883104,80.4834737271399]
};		

t_best_OUFIM={[7.63437186283769,15.8005540434218,33.3382250451239]
[6.61404564863434e-07,8.36271649838249,16.5826941192415,78.7187814732675]%34.7187814732675
[7.82033147903402e-07,8.24621981659480,16.4095729804949,30.0111305051334,79.9999965854638]
[6.01187288621078e-07,7.22222284904399,12.9516044663223,18.3837728715959,31.2295955995076,79.9999967941538]
[5.22239927753310e-07,6.16915212420954,10.4399825971812,15.2717163062060,20.2232787410545,32.1979680046146,79.9999968668265]
[4.98576739422388e-07,5.84341152554804,9.77027761173052,14.2503755818218,18.4214379782159,25.6610483427347,36.2301652737999,79.9999969978647]
[1.90484082558543e-06,5.29293551495518,8.75771126292488,12.4152074349164,15.9969137708275,19.9290581458478,27.4532981863932,38.1092767687901,79.9999878481224]
[1.85868261715652e-06,4.83231839740851,7.97708652282505,11.0480101556969,14.3638164678974,17.5188346247261,21.6149491589979,28.9236100909354,39.9099730205308,79.9999885830239]
};  

set_numpoints=[3,4,5,6,7,8,9,10];
ave_num=1000;
SigmaOU=0.16;%%%%consistent with OU noise
phi=0.02;
r=0.2; K=50; N0=4.5;
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'C_0'};   
l_set=length(t_best_OUFIM); 
t_random=cell(l_set,1);
t_even=cell(l_set,1);
 
for i=1:l_set
    n_s=set_numpoints(i); 
    t_random{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));
    t_even{i}=linspace(t_initial, TF, n_s);%equally spaced points
end

% t_vec={t_random,t_even,t_best_IIDFIM};
t_vec={t_even,t_best_IIDFIM,t_best_OUFIM};
% prm_interest={'N_0','K'}; 
% prm_interest={'N_0'}; 
% prm_interest={'K'}; 
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

%%%%%%%%%%%parameters for OU process noise
 Sigma=sqrt(SigmaOU^2/(2*phi));
% lb_p=[0.1,25,0.1]; % the bounds of the parameters r, K and N0 for computing profile likehood.
% ub_p=[0.4,75,15];
 lb_p=[0.1,40,1]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.3,60,10];   

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
             measurements=N_real+Generate_exact_OU(phi,SigmaOU,t_opt); % OU measurement noise
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
%plot the figures
 CI_bounds=[0.015,1,1;
                        0.065,3.7,3.5];

 fig=figure('Position',[20 20 1000 320],'color','w');
  %   [left bottom width height]
for j=1:num_param
%       subplot('Position',[0.06+(j-1)*0.31,0.12,0.265,0.75]);
       subplot('Position',[0.06+(j-1)*0.31,0.13,0.265,0.75]);
%     subplot(1,num_param,j,'Position',[0.05+(j-1)*0.3,0.12,0.26,0.75]);
%       subplot(1,num_param,j);
%         plot(set_numpoints, CI_w(kk,:,j+3), 'o--', set_numpoints, CI_w(kk,:,j), 'd--',  'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);
    if j==1
        plot(set_numpoints(3:end),  CI_w(3:end,j+0),'*-',set_numpoints,CI_w(:,j+num_param*(num_types-2)), 'x:',...
                 set_numpoints, CI_w(:,j+num_param*(num_types-1)), 'o--', ...
                 'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);
    else
        plot(set_numpoints,  CI_w(:,j+0),'*-',set_numpoints,CI_w(:,j+num_param*(num_types-2)), 'x:',...
            set_numpoints, CI_w(:,j+num_param*(num_types-1)), 'o--', ...
            'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);
    end
    xlim([set_numpoints(1),set_numpoints(end)]);
     ylim([CI_bounds(1,j),CI_bounds(2,j)]);
    xlabel('Number of Points');
    ylabel('Confifence Intervals Width');
       legend('IID FIM','OU FIM','OU Sobol','Location','Northeast');
%   legend('Evenly','Optimized FIM','Optimized Sobol'); 
   
    title(sprintf('Parameter %s', prm_name{j}));   
end



   
elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
