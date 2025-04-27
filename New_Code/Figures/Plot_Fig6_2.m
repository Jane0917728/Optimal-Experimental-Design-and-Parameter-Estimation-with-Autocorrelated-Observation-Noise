clear all;
close all;

tic

t_best_OUFIM={[7.63437186283769,15.8005540434218,33.3382250451239]
    [6.61404564863434e-07,8.36271649838249,16.5826941192415,78.7187814732675]%34.7187814732675
    [7.82033147903402e-07,8.24621981659480,16.4095729804949,30.0111305051334,79.9999965854638]
    [6.01187288621078e-07,7.22222284904399,12.9516044663223,18.3837728715959,31.2295955995076,79.9999967941538]
    [5.22239927753310e-07,6.16915212420954,10.4399825971812,15.2717163062060,20.2232787410545,32.1979680046146,79.9999968668265]
    [4.98576739422388e-07,5.84341152554804,9.77027761173052,14.2503755818218,18.4214379782159,25.6610483427347,36.2301652737999,79.9999969978647]
    [1.90484082558543e-06,5.29293551495518,8.75771126292488,12.4152074349164,15.9969137708275,19.9290581458478,27.4532981863932,38.1092767687901,79.9999878481224]
    [1.85868261715652e-06,4.83231839740851,7.97708652282505,11.0480101556969,14.3638164678974,17.5188346247261,21.6149491589979,28.9236100909354,39.9099730205308,79.9999885830239]
    };

t_best_OUSob={[10,16,30]
    [2,10,16,30]
    [2,10,16,28,80]
    [0,8,12,16,28,80]
    [0,8,12,16,20,30,80]
    [0,8,10,14,18,26,36,80]
    [0,6,10,14,16,20,26,36,80]
    [0,4,8,10,14,16,20,26,36,80]};

set_numpoints=[3,4,5,6,7,8,9,10];
ave_num=1000;
Sigma=0.16; % onsistent with OU noise
r=0.2;
K=50;
N0=4.5;
phi=0.02;
t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;
prm_all=[r K N0]'; 
prm_name_all = {'r', 'K', 'C_0'};
l_set=length(t_best_OUFIM);
t_random=cell(l_set,1);
t_even=cell(l_set,1);

for i=1:l_set
    n_s=set_numpoints(i);
    t_random{i}=sort(t_initial+(TF-t_initial)*rand(1,n_s));
    t_even{i}=linspace(t_initial, TF, n_s); % equally spaced points
end

t_vec={t_even,t_best_OUFIM,t_best_OUSob};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

% parameters for OU process noise
lb_p=[0.1,40,1]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.3,60,10];
numpts=70; % parameter number of the profile likelihood
num_types=length(t_vec); % number of types, such as OU noise with FIM, or OU noise with Sob index.
CI_w=zeros(l_set,num_param*num_types);%tThe first  means the number of change parameters
CI=zeros(l_set*4,num_param*num_types); % 4 denotes all information CI_lower bound, estimated value,
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

    disp(j)

    t_v=t_vec{j};
    for k=1:l_set
        t_opt=t_v{k,1};
        N_real=Nfunction(r,K,t_opt,N0);
        n_s=length(t_opt);
        SIG=zeros(n_s,n_s); % The covariance matrix
        for s = 1:n_s
            for tt = 1:n_s
                SIG(s, tt) = Sigma^2*exp(-phi*abs(t_opt(tt)-t_opt(s)))/(2*phi);
            end
        end

        parfor i=1:ave_num
            measurements=N_real+Generate_exact_OU(phi,Sigma,t_opt);
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
        % average only for the results of
        % successfully getting the confident intervel.
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

save('Plot_Fig6_2.mat');

elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);

%% 

clear all;
close all;

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

load('Plot_Fig6_2mat')

fig=figure('Position',[20 20 1000 320],'color','w');
for j=1:num_param
    subplot('Position',[0.06+(j-1)*0.31,0.14,0.265,0.75]);
    if j==1
        plot(set_numpoints(3:end),CI_w(3:end,j+0),'o-',...
            set_numpoints,CI_w(:,j+num_param*(num_types-2)),'o--',...
            set_numpoints,CI_w(:,j+num_param*(num_types-1)),'o-.', ...
            'MarkerSize',6,'LineWidth',2.0);
    else
        plot(set_numpoints,CI_w(:,j+0),'o-',...
            set_numpoints,CI_w(:,j+num_param*(num_types-2)),'o--',...
            set_numpoints,CI_w(:,j+num_param*(num_types-1)),'o-.', ...
            'MarkerSize',6,'LineWidth',2.0);
    end
    if j==1
        ytickformat('%.2f');
        ylim([0.0 0.06])
        yticks([0.00 0.02 0.04 0.06]);
        title('$r$','Interpreter', 'latex')
        ylabel('Mean confidence interval width');
    end
    if j==2
        ytickformat('%.1f');
        ylim([0 4])
        title('$K$','Interpreter', 'latex')
        yticks([0.0 1.0 2.0 3.0 4.0]);
    end
    if j==3
        ytickformat('%.1f');
        ylim([0 4])
        title('$C_0$','Interpreter', 'latex')
        yticks([0.0 1.0 2.0 3.0 4.0]);
    end
    xlabel('Number of points');
    if  j==2
        legend('Evenly','Optimized FIM','Optimized Sobol');
    end
end

print('Fig6_2','-dpng','-r300')