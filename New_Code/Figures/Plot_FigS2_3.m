clear all
close all

t_best_IIDSob={[8,16,80]
    [8,16,18,80]
    [8,10,16,18,80]
    [8,10,16,18,74,80]
    [8,10,16,18,74,76,80]
    [8,10,14,16,18,74,76,80]
    [8,10,12,16,18,20,74,76,80]
    [8,10,12,16,18,20,74,76,78,80]};

t_best_IIDFIM={[6.34994863705261,16.7861782858528,80.4979392777307]
    [6.20977450676229,15.2750903927458,18.0692952336378,80.4922051169937]
    [6.27577407941707,15.7021317824451,17.7021364588466,78.4963945183739,80.4986631701891]
    [5.34126363308649,7.34126604902159,15.8463055940299,17.8463085938521,78.4948727676437,80.4962916056809]
    [5.33730820291159,7.33732568473259,15.8329599735807,17.8329799902337,76.4053885595053,78.4511924475796,80.4786168679362]
    [5.19742621425561,7.19742714665278,14.7063616024953,16.7063635963056,18.7063650782529,76.4943377507751,78.4971540817235,80.4986760854827]
    [4.34607220991320,6.34608563188164,8.34610109975754,14.9428091241092,16.9428267679213,18.9428441988758,76.4077382362805,78.4524124723076,80.4791442100305]
    [4.33685461950222,6.33686782322273,8.33688291127786,14.9132214464828,16.9132380161570,18.9132539057305,74.3998089588690,76.4404162632139,78.4644728883104,80.4834737271399]
    };

t_best_OUSob={[8,16,72]
[8,16,56,80]
[8,16,46,64,80]
[8,14,18,44,62,80]
[8,14,18,40,52,66,80]
[8,12,16,20,38,52,68,80]
[8,12,16,20,34,44,56,68,80]
[8,12,16,20,32,42,52,60,70,80]};
    
t_best_OUFIM={[6.36724130734275,16.7688695158861,79.9982905015772]
[6.36397256137513,16.7639805367088,56.9622958881232,79.9999551302772]
[5.18594144372386,12.4114706755654,18.9437869989987,57.2117104889478,79.9999464965157]
[5.06374516742977,12.0852900761270,18.6470226931543,46.6586677915456,63.3981168384554,79.9999882077699]
[3.41350476819170,8.75718073974524,15.0060117923993,20.5951799903851,46.9429332072944,63.5394376792389,79.9999888046925]
[3.35004450071383,8.65050028345085,14.8484100616149,20.2982541688903,40.4721542134551,53.8333827432587,66.9283423404579,79.9999946994061]
[2.47475433649929,7.38275234443300,12.5097871483522,17.2215957631767,22.6999703809768,41.1504007580529,54.2704617423436,67.1462847267834,79.9999947370060]
[2.30146422114702,7.14704346032010,12.0937164371786,16.7721720467986,21.9071438260459,36.3160435162458,47.5682316619041,58.4173930674881,69.2123845913418,79.9999965470977]
  };  

set_numpoints=[3,4,5,6,7,8,9,10];
ave_num=1000;

r=0.2;
K=50;
N0=4.5;
phi=0.3;
SigmaC=2;
Sigma=sqrt(SigmaC^2*(2*phi));

t_interval=2;
t_initial=0;
TF=80+t_initial;
t=t_initial:t_interval:floor(TF/t_interval)*t_interval+t_initial;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'C_0'};
l_set=length(set_numpoints); % number of set_numberpoints
t_vec={t_best_IIDFIM,t_best_IIDSob,t_best_OUFIM,t_best_OUSob};

prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);
r_index=1;
K_index=1;
N0_index=1;

lb_p=[0.05,20,0.1];
ub_p=[0.4,70,30];
numpts=70; % parameter number of the profile likelihood

num_types=length(t_vec); % number of types, such as OU noise with FIM, or OU noise with Sob index.
CI_w=zeros(l_set,num_param*num_types); % num_param*num_types means number of parameters times number
% number of types
CI=zeros(l_set*4,num_param*num_types); % 4 denotes all information CI_lower bound, estimated value,
% CI upper bound and width of CI.
CI_1=CI;
CI_w1=CI_w;
CI_width=zeros(l_set,num_param*num_types);

CI_param=CI_width;

suc=zeros(l_set,num_param*num_types); % success rate among the total ave_num runs.

CI_r=zeros(3,ave_num); % CI_r denotes the lower bound, the estimated value and the upper bound.
CI_K=CI_r; % CI_K denotes the lower bound, the estimated value and the upper bound.
CI_N0=CI_r; % CI_N0 denotes the lower bound, the estimated value and the upper bound.

flag=false;
for j=1:num_types
    j
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
        SIG_OU=SIG;
        SIG_IID=eye(size(SIG))*SigmaC^2;

        parfor i=1:ave_num
            measurements=N_real+Generate_exact_OU(phi,Sigma,t_opt);
            [CI_r(:,i),CI_K(:,i),CI_N0(:,i),MSE_est(i)] = ...
                get_CI(r,K,N0,numpts,SIG_OU,lb_p,ub_p,measurements,t_opt,num_param);
        end

        %%
        % The average of all confidence intervals, including the unsuccessful ones.
        CI_r_ave1 = mean(CI_r,2);
        CI_K_ave1= mean(CI_K,2);
        CI_N0_ave1 = mean(CI_N0,2);

        % Indentify the indices corresponding to valid confidence
        r_index=(CI_r(1,:)>-0.2)&(CI_r(3,:)<1);
        K_index=(CI_K(1,:)>10)&(CI_K(3,:)<110);
        N0_index=(CI_N0(1,:)>-1)&(CI_N0(3,:)<20);

        % Check r_index and compute mean if not empty
        if any(r_index)
            CI_r_ave = mean(CI_r(:, r_index), 2);
        else
            CI_r_ave = ones(size(CI_r, 1), 1)*1.2; % Or set to another appropriate value
            disp('r All fail');
            disp(CI_r);
        end

        % Check K_index and compute mean if not empty
        if any(K_index)
            CI_K_ave = mean(CI_K(:, K_index), 2);
        else
            CI_K_ave = ones(size(CI_K, 1), 1)*100; % Or set to another appropriate value
            disp('K All fail');
            disp(CI_K);
        end

        % Check N0_index and compute mean if not empty
        if any(N0_index)
            CI_N0_ave = mean(CI_N0(:, N0_index), 2);
        else
            CI_N0_ave = ones(size(CI_N0, 1), 1)*21; % Or set to another appropriate value
            disp('C_0 All fail');
            disp(CI_N0);
        end

        % average only for the results of
        CI_f_r1=[CI_r_ave1;CI_r_ave1(3)-CI_r_ave1(1)];
        CI_f_K1=[CI_K_ave1;CI_K_ave1(3)-CI_K_ave1(1)];
        CI_f_N01=[CI_N0_ave1;CI_N0_ave1(3)-CI_N0_ave1(1)];
        CI_f_r=[CI_r_ave;CI_r_ave(3)-CI_r_ave(1)];
        CI_f_K=[CI_K_ave;CI_K_ave(3)-CI_K_ave(1)];
        CI_f_N0=[CI_N0_ave;CI_N0_ave(3)-CI_N0_ave(1)];

        CI_1(4*(k-1)+1:4*k,((j-1)*num_param+1):j*num_param)=[CI_f_r1,CI_f_K1,CI_f_N01]; % CI
        % CI_width and CI_param, where k=1:l_set number of typeofpoints,
        % j=1:number of types
        CI_w1(k,((j-1)*num_param+1):j*num_param)=CI_1(4*k,((j-1)*num_param+1):j*num_param);
        CI(4*(k-1)+1:4*k,((j-1)*num_param+1):j*num_param)=[CI_f_r,CI_f_K,CI_f_N0];%CI
        % CI_width and CI_param
        CI_w(k,((j-1)*num_param+1):j*num_param)=CI(4*k,((j-1)*num_param+1):j*num_param);

        CI_param(k,((j-1)*num_param+1):j*num_param)=CI(4*(k-1)+2,((j-1)*num_param+1):j*num_param);
        % compute the Identification success ratio
        suc(k,((j-1)*num_param+1):j*num_param)=[sum(r_index)/ave_num,sum(K_index)/ave_num,sum(N0_index)/ave_num];

    end

end

CI_all=CI;
CI_parameters=CI_param;
suc_times=suc*ave_num;

save('Plot_FigS2_3.mat');

elapsed_time = toc;
disp(['Elapsed time: ',num2str(elapsed_time),' seconds']);

%% 

clear all
close all

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

load('Plot_FigS2_3.mat');

fig=figure('Position',[20 20 1000 320],'color','w');
index=1:8;
index2=[4,8,12,16,20,24,28,32];
CI_w=CI(index2,:);
for j=1:num_param
    subplot('Position',[0.06+(j-1)*0.32,0.14,0.265,0.75]);
    plot(set_numpoints(index),CI_w(index,j),'o-',...
        set_numpoints(index),CI_w(index,j+num_param*1),'o--',...
        set_numpoints(index),CI_w(index,j+num_param*2),'o-.', ...
        set_numpoints(index), CI_w(index,j+num_param*3), 'o:', ...
        'MarkerSize',6,'LineWidth',2);
    xlim([set_numpoints(index(1)),set_numpoints(index(end))]);
    xlabel('Number of points');
    if j==1
        ytickformat('%.2f');
        ylim([0 0.15])
        yticks([0.00 0.05 0.10 0.15]);
        title('$r$','Interpreter', 'latex')
        ylabel('Mean confidence interval width');
    end
    if j==2
        ytickformat('%.0f');
        ylim([0 10])
        yticks([0 5 10]);
        title('$K$','Interpreter', 'latex')
        legend('IID FIM','IID Sobol','OU FIM','OU Sobol','Location','SouthWest');
    end

    if j==3
        ytickformat('%.0f');
        ylim([0 8])
        yticks([0 4 8]);
        title('$C_0$','Interpreter', 'latex')
    end
end

print('FigS2_3','-dpng','-r300')

