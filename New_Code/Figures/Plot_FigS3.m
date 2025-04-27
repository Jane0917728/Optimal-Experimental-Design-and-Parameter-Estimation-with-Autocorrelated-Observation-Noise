clear all
close all

tic

set_numpoints=[5,6,7,8,9];
ave_num=1000;
Sigma=0.8; % consistent with OU noise Simga_OU=0.6 as phi=0.02
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
    t_even{i}=linspace(t_initial, TF, n_s); % equally spaced points
end

t_check={[0,20,40,60,80]
    [0,20,40,60,64,80]
    [0,13.3,20,40,60,80]
    [0,11.4,20,40,45.7,60,68.6,80]
    [0,10,20,30,40,50,60,70,80]
    };

t_vec={t_even,t_check};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

% parameters for OU process noise
lb_p=[0.05,25,0.1]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.4,75,15];

CI_bounds=[0.015,1,1;
    0.065,3.7,3.5];

numpts=70; % parameter number of the profile likelihood
num_types=length(t_vec); % number of types, such as OU noise with FIM, or OU noise with Sob index.
CI_w=zeros(l_set,num_param*num_types); % the first  means the number of change parameters
CI=zeros(l_set*4,num_param*num_types); % 4 denotes all information CI_lower bound, estimated value,
% CI upper bound and width of CI.

CI_width=zeros(l_set,num_param*num_types);
CI_param=CI_width;
MSE_param=CI_w; % MSE of each parameters
suc=zeros(l_set,num_param*num_types); % 8*4
MSE_est_ave=zeros(l_set,num_param*num_types); % the MSE of the estimated parameters
CI_r=zeros(3,ave_num); % CI_r denotes the lower bound, the estimated value and the upper bound.
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
            measurements=N_real+ Sigma*randn(size(N_real)); % IID measurement noise
            [CI_r(:,i),CI_K(:,i),CI_N0(:,i),MSE_est(i)] = ...
                get_CI(r,K,N0,numpts,SIG,lb_p,ub_p,measurements,t_opt,num_param);
        end

        r_index=(CI_r(1,:)>-0.2)&(CI_r(3,:)<1);
        K_index=(CI_K(1,:)>10)&(CI_K(3,:)<110);
        N0_index=(CI_N0(1,:)>-1)&(CI_N0(3,:)<20);
        % average for all the results of ave_num
        % running
        CI_r_ave1 = mean(CI_r,2);
        CI_K_ave1= mean(CI_K,2);
        CI_N0_ave1 = mean(CI_N0,2);

        CI_r_ave = mean(CI_r(:,r_index),2);
        CI_K_ave= mean(CI_K(:,K_index),2);
        CI_N0_ave = mean(CI_N0(:,N0_index),2);
        
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

save('Plot_FigS3.mat');

elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);

%% 

clear all
close all

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

load('Plot_FigS3.mat')

fig=figure('Position',[20 20 1000 320],'color','w');
for j=1:num_param
    subplot('Position',[0.06+(j-1)*0.31,0.14,0.265,0.75]);
    plot(set_numpoints(1:end),CI_w(1:end,j+0),'o-',...
        set_numpoints,CI_w(:,j+num_param*(num_types-1)),'o--',...
        'MarkerSize',6,'LineWidth',2.0);
    xlim([set_numpoints(1),set_numpoints(end)]);
    
    if j==1
        ytickformat('%.2f');
        ylim([0.0 0.08])
        yticks([0.00 0.02 0.04 0.06 0.08]);
        title('$r$','Interpreter', 'latex')
        ylabel('Mean confidence interval width');
    end
    if j==2
        ytickformat('%.1f');
        ylim([0 2])
        ytickformat('%.1f');
        title('$K$','Interpreter', 'latex')
        yticks([0.0 0.5 1.0 1.5 2.0]);
        legend('Evenly','Test');
    end
    if j==3
        ytickformat('%.0f');
        ylim([0 4])
        title('$C_0$','Interpreter', 'latex')
        yticks([0.0 1.0 2.0 3.0 4.0]);
    end
    xlabel('Number of points');
    if  j==2
        legend('Evenly','Test');
    end
    xlabel('Number of points');
end
print('FigS3','-dpng','-r300')
