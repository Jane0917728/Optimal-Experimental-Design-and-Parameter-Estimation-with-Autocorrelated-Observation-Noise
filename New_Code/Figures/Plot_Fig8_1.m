clear all;
close all;

tic

load('Plot_Fig8_FIM_t_best.mat');
load('Plot_Fig8_Sob_t_best.mat');

phi=[0.0100
    0.0206
    0.0424
    0.0874
    0.1800
    0.2300
    0.5250
    0.8200
    1.1150
    1.4100
    1.7050
    2.0000];

% t_best_OUFIM={[1.08436537138667e-06,4.79779544093355,7.83949267092782,10.7637071177565,14.0616585716086,17.1200694942515,20.8420081309374,27.5609838074098,36.3907380846756,79.9999770533238]
% [1.87374382812792e-06,4.83231345613365,7.97714665452729,11.0481364715128,14.3637854910520,17.5187698386371,21.6149227863610,28.9235154074201,39.9100252924649,79.9999886838623]
% [1.23476166751385e-05,5.25856267301390,8.82142429998489,12.6084315029357,16.1647967749064,20.2194176258275,28.8359979292716,40.9702983300488,60.0962814686171,79.9999917120394]
% [3.13931079786280,7.35797228620882,11.3150926748034,15.2930829663892,19.3903524915685,28.1957131351803,39.0305946751084,52.3628037213324,66.1598936207319,79.9999982142682]
% [3.06672148791616,7.47507007176796,11.7826760478282,15.9080719560237,20.4015027479486,31.6339291648082,43.3041479021851,55.4769131064347,67.7334182501756,79.9999985502799]
% [2.75341295250315,7.37408277146741,11.9083042562720,16.2085970684083,20.9983891476464,32.9346198646150,44.5955203694181,56.3777986216219,68.1870104699910,79.9999987321380]
% [2.39660247569670,7.24696227481680,12.0523829288142,16.5806442521275,21.7126645289012,34.2776642877302,45.8296156193570,57.2301104531145,68.6158001293733,79.9999980324751]
% [2.52525686927827,7.13216632240105,11.9859297110471,16.6186241906652,21.4187236516069,39.4366311561430,49.9163785486568,59.9886559647637,69.9991686249221,79.9999955935505]
% [2.55491251801156,6.16166385327941,9.76661338347361,14.0307830515877,17.5984032318381,21.4503660225014,51.8685311355067,61.4145071240480,70.7264322606173,79.9998759363814]
% [3.21775903206623,6.31196399818162,9.52997944910866,13.9718490731083,17.1750279014725,20.4143359668692,56.6519897828257,64.5886884735500,72.3162892433386,79.9997350260320]
% [3.84021409435062,6.32790469534669,8.91557295009929,14.3708526952945,16.9839637650038,19.5697185940709,61.7849487133504,67.9970257405257,74.0238176533114,79.9998388492052]
% [4.22826820588485,6.31794836813069,8.46600639397794,14.7464184704408,16.9286382515405,19.0935510780972,65.0447506275423,70.1570645506861,75.1047759195942,79.9997779763455]};

% t_best_OUSob={[0,4,8,10,14,16,18,24,32,80]
% [0,4,8,10,14,16,20,26,34,80]
% [0,8,10,14,16,20,26,36,58,80]
% [0,8,10,14,18,26,36,50,66,80]
% [8,10,14,18,28,38,50,60,70,80]
% [8,12,16,20,32,40,50,60,70,80]
% [6,10,14,16,20,40,48,58,70,80]
% [6,8,10,14,16,20,52,60,70,80]
% [6,8,10,14,16,18,56,64,72,80]
% [6,8,10,14,16,18,58,64,72,80]
% [6,8,10,14,16,18,62,68,74,80]
% [6,8,10,14,16,18,62,68,74,80]
% };

ave_num=1000;
r=0.2; 
K=50; 
N0=4.5;
t_initial=0;
TF=80+t_initial;

prm_all=[r K N0]'; 
prm_name_all = {'r', 'K', 'C_0'};

SigmaC=2;
num_phi=length(phi);
l_set=num_phi;
t_vec={t_best,t_bestSob};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

lb_p=[0.05,25,0.5]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.65,75,14];

numpts=70; % parameter number of the profile likelihood
num_types=length(t_vec); % number of types, such as OU noise with FIM, or OU noise with Sob index.
CI_w=zeros(l_set,num_types*num_param); % the first  means the number of change parameters
% the average without the practical unidentifiable
CI=zeros(l_set*4,num_types*num_param); % 4 denotes all information CI_lower bound, estimated value,
% CI upper bound and width of CI.
CI_1=CI;
CI_w1=CI_w;
CI_width=zeros(l_set,num_types*num_param);
CI_param=CI_width;
suc=zeros(l_set,num_types*num_param); %8*4
CI_r=zeros(3,ave_num);  % CI_r denotes the lower bound, the estimated value and the upper bound.
CI_K=CI_r; % CI_K denotes the lower bound, the estimated value and the upper bound.
CI_N0=CI_r; % CI_N0 denotes the lower bound, the estimated value and the upper bound.

Sigma=sqrt(SigmaC^2*(2*phi));
for k=1:num_types

    disp(k)
    
    t_opt1=t_vec{k};
    
    for j=1:l_set
        
        t_opt=t_opt1{j};
        N_real=Nfunction(r,K,t_opt,N0);
        n_s=length(t_opt);
        SIG=zeros(n_s,n_s); % The covariance matrix

        for s = 1:n_s
            for tt = 1:n_s
                SIG(s, tt) = Sigma(j)^2/(2*phi(j))*exp(-phi(j)*abs(t_opt(tt)-t_opt(s)));
            end
        end
        SIG_OU=SIG;
        phi1=phi(j);
        Sigma1=Sigma(j);
        parfor i=1:ave_num
            measurements=N_real+Generate_exact_OU(phi1,Sigma1,t_opt);
            [CI_r(:,i),CI_K(:,i),CI_N0(:,i),MSE_est(i)] = ...
                get_CI(r,K,N0,numpts,SIG_OU,lb_p,ub_p,measurements,t_opt,num_param);
        end
        r_index=(CI_r(1,:)>-0.2)&(CI_r(3,:)<1);
        K_index=(CI_K(1,:)>10)&(CI_K(3,:)<110);
        N0_index=(CI_N0(1,:)>-1)&(CI_N0(3,:)<20);
        
        % average for all the results of ave_num
        % running
        CI_r_ave1 = mean(CI_r,2);
        CI_K_ave1= mean(CI_K,2);
        CI_N0_ave1 = mean(CI_N0,2);

        % Check r_index and compute mean if not empty
        if any(r_index)
            CI_r_ave = mean(CI_r(:, r_index), 2);
        else
            CI_r_ave = ones(size(CI_r, 1), 1)*1.2; % Or set to another appropriate value
            disp('All fail');
            disp(CI_r);
        end

        % Check K_index and compute mean if not empty
        if any(K_index)
            CI_K_ave = mean(CI_K(:, K_index), 2);
        else
            CI_K_ave = ones(size(CI_K, 1), 1)*100; % Or set to another appropriate value
            disp('All fail');
            disp(CI_K);
        end

        % Check N0_index and compute mean if not empty
        if any(N0_index)
            CI_N0_ave = mean(CI_N0(:, N0_index), 2);
        else
            CI_N0_ave = ones(size(CI_N0, 1), 1)*21; % Or set to another appropriate value
            disp('All fail');
            disp(CI_K);
        end

        CI_f_r1=[CI_r_ave1;CI_r_ave1(3)-CI_r_ave1(1)];
        CI_f_K1=[CI_K_ave1;CI_K_ave1(3)-CI_K_ave1(1)];
        CI_f_N01=[CI_N0_ave1;CI_N0_ave1(3)-CI_N0_ave1(1)];
        CI_f_r=[CI_r_ave;CI_r_ave(3)-CI_r_ave(1)];
        CI_f_K=[CI_K_ave;CI_K_ave(3)-CI_K_ave(1)];
        CI_f_N0=[CI_N0_ave;CI_N0_ave(3)-CI_N0_ave(1)];

        CI_1(4*(j-1)+1:4*j,((k-1)*num_param+1):k*num_param)=[CI_f_r1,CI_f_K1,CI_f_N01];%CI
        % CI_width and CI_param
        CI_w1(j,((k-1)*num_param+1):k*num_param)=CI_1(4*j,((k-1)*num_param+1):k*num_param);
        CI(4*(j-1)+1:4*j,((k-1)*num_param+1):k*num_param)=[CI_f_r,CI_f_K,CI_f_N0];%CI
        % CI_width and CI_param
        CI_w(j,((k-1)*num_param+1):k*num_param)=CI(4*j,((k-1)*num_param+1):k*num_param);

        CI_param(j,((k-1)*num_param+1):k*num_param)=CI(4*(j-1)+2,((k-1)*num_param+1):k*num_param);
        % compute the Identification success ratio
        suc(j,((k-1)*num_param+1):k*num_param)=[sum(r_index)/ave_num,sum(K_index)/ave_num,sum(N0_index)/ave_num];
    end
end

save('Plot_Fig8_1.mat')

elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);

%% 

clear all;
close all;

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

load('Plot_Fig8_1.mat')

fig=figure('Position',[20 20 1000 320],'color','w');

for k=1:num_param
    subplot('Position',[0.06+(k-1)*0.31,0.14,0.265,0.75]);
    plot(phi,CI_w(:,k),'o-', phi,CI_w(:,k+num_param*1),'o--','MarkerSize',6,'LineWidth',2.0);
    xlabel('$\phi$', 'Interpreter', 'latex');
    if k==1
        ylabel('Mean confidence interval width');
        title('$r$','Interpreter', 'latex')
        xticks([0 0.5 1.0 1.5 2.0]);
        xtickformat('%.1f');
        ytickformat('%.2f');
        ylim([0 0.15])
        yticks([0.00 0.05 0.1 0.15]);
    end
    if k==2
        title('$K$','Interpreter', 'latex')
        legend('OU FIM','OU Sobol');
        xticks([0 0.5 1.0 1.5 2.0]);
        xtickformat('%.1f');
        ylim([0 8])
        yticks([0 2 4 6 8]);
    end
    if k==3
        ylim([0 6])
        yticks([0 2 4 6])
        title('$C_0$','Interpreter', 'latex')
        xticks([0 0.5 1.0 1.5 2.0]);
        xtickformat('%.1f');
    end
end
print('Fig8_1','-dpng','-r300')
