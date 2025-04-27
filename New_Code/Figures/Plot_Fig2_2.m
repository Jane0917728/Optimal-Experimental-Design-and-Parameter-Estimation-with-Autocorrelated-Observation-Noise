clear all;
close all;

set(0,'DefaultAxesFontSize',13);
set(0,'DefaultTextFontSize',13);

ave_num=1000;
r=0.2; 
K=50; 
N0=4.5;
t_interval=2;
t_initial=0;
TF=80+t_initial;
prm_all=[r K N0]'; 
prm_name_all = {'r', 'K', 'C_0'};
n_s=11; % number of measurement points
t=linspace(t_initial,TF,n_s);
num_sigma=8;
phi=0.1;

Sigma=linspace(0.5,1.5,num_sigma)';
sigmaOU=Sigma;
SigmaEqual=sqrt(sigmaOU.^2/(2*phi));
sigmaIID=Sigma;
l_set=num_sigma;

prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

% parameters for OU process noise

lb_p=[0.05,25,0.5]; % the bounds of the parameters r, K and N0 for computing profile likehood.
ub_p=[0.65,75,14];
numpts=60; % parameter number of the profile likelihood
num_types=3; % number of types, such as OU noise, IID noise and misspecified noise.
CI_w=zeros(l_set,num_types*num_param); % the first means the number of change parameters
%the average without the practical unidentifiable
CI=zeros(l_set*4,num_types*num_param); % 4 denotes all information CI_lower bound, estimated value,
%CI upper bound and width of CI.
CI_1=CI;
CI_w1=CI_w;
CI_width=zeros(l_set,num_types*num_param);

CI_param=CI_width;
MSE_param=CI_w; % MSE of each parameters
suc=zeros(l_set,num_types*num_param); %8*4
MSE_est_ave=zeros(l_set,num_types*num_param); % the MSE of the estimated parameters
CI_r=zeros(3,ave_num);  % CI_r denotes the lower bound, the estimated value and the upper bound.
CI_K=CI_r; % CI_K denotes the lower bound, the estimated value and the upper bound.
CI_N0=CI_r; % CI_N0 denotes the lower bound, the estimated value and the upper bound.
MSE_est=zeros(ave_num,1); % record the mean of square error between the estimated value and the real value.

OU_noise=zeros(ave_num,n_s);
t_opt=t;
N_real=Nfunction(r,K,t_opt,N0);
SIG_OU=zeros(n_s,n_s); % The covariance matrix
for j=1:l_set

    disp(j)

    for s = 1:n_s
        for tt = 1:n_s
            SIG_OU(s, tt) = sigmaOU(j)^2/(2*phi)*exp(-phi*abs(t_opt(tt)-t_opt(s)));
        end
    end
    SIG_IID=eye(n_s,n_s)*sigmaIID(j)^2;

    for na=1:ave_num
        OU_noise(na,:)=Generate_exact_OU(phi,sigmaOU(j),t_opt);
    end
    IID_noise=sigmaIID(j)*randn(ave_num,n_s);

    for k=1:num_types
        if k==1
            %OU noise
            SIG=SIG_OU;
            noise=OU_noise;
        end
        if k==2
            %IID noise
            SIG=SIG_IID;
            noise=IID_noise;
        end
        if k==3
            %Misspecified noise
            SIG=SIG_IID;
            noise=OU_noise;
        end
        parfor i=1:ave_num
            measurements=N_real+noise(i,:);
            [CI_r(:,i),CI_K(:,i),CI_N0(:,i),MSE_est(i)] = ...
                get_CI(r,K,N0,numpts,SIG,lb_p,ub_p,measurements,t_opt,num_param);
          
        end
        r_index=(CI_r(1,:)>-0.2)&(CI_r(3,:)<1);
        K_index=(CI_K(1,:)>10)&(CI_K(3,:)<110);
        N0_index=(CI_N0(1,:)>-1)&(CI_N0(3,:)<20);
        
        CI_r_ave1 = mean(CI_r,2);
        CI_K_ave1= mean(CI_K,2);
        CI_N0_ave1 = mean(CI_N0,2);
        
        % Check r_index and compute mean if not empty
        if any(r_index)
            CI_r_ave = mean(CI_r(:, r_index), 2);
        else
            CI_r_ave = ones(size(CI_r, 1), 1)*1.2; % Or set to another appropriate value
            disp('r all fail');
            disp(CI_r);
        end

        % Check K_index and compute mean if not empty

        if any(K_index)
            CI_K_ave = mean(CI_K(:, K_index), 2);
        else
            CI_K_ave = ones(size(CI_K, 1), 1)*100; % Or set to another appropriate value
            disp('K all fail');
            disp(CI_K);
        end

        % Check N0_index and compute mean if not empty
        if any(N0_index)
            CI_N0_ave = mean(CI_N0(:, N0_index), 2);
        else
            CI_N0_ave = ones(size(CI_N0, 1), 1)*21; % Or set to another appropriate value
            disp('C_0 all fail');
            disp(CI_N0);
        end
        
        MSE_ave_r=sqrt(mean((CI_r(2,:)-r).^2,2));
        MSE_ave_K=sqrt(mean((CI_K(2,:)-K).^2,2));
        MSE_ave_N0=sqrt(mean((CI_N0(2,:)-N0).^2,2));
        
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
        MSE_param(j,((k-1)*num_param+1):k*num_param)=[MSE_ave_r,MSE_ave_K,MSE_ave_N0];
        suc(j,((k-1)*num_param+1):k*num_param)=[sum(r_index)/ave_num,sum(K_index)/ave_num,sum(N0_index)/ave_num];
        
    end
end

%%

save('Plot_Fig2_2mat');

%%

clear all
close all

load('Plot_Fig2_2.mat')

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

figure(1);
hFig = gcf;
set(hFig,'Position',[20 20 1000 320],'color','w');

for k=1:num_param
    subplot('Position',[0.06+(k-1)*0.31,0.14,0.265,0.75]);
    plot(Sigma,CI_w(:,k),'-o',Sigma,CI_w(:,k+num_param*1),'--o', Sigma, CI_w(:,k+num_param*2),'-.o','MarkerSize',6,'LineWidth',2);
    xlim([Sigma(1),Sigma(end)]);
    xlabel('\sigma');
    if k==1
        ylabel('Mean confidence interval width');
    end
    if k==1
        ylim([0 0.25])
        xtickformat('%.1f');
        ytickformat('%.2f');
        xticks([1.0 2.0 3.0 4.0]);
        yticks([0.00 0.05 0.10 0.15 0.20 0.25]);
        title('$r$','Interpreter', 'latex')
    end
    if k==2
        ylim([0 25])
        xtickformat('%.1f');
        ytickformat('%.0f');
        xticks([1.0 2.0 3.0 4.0]);
        % yticks([0 5 10 15 20 25]);
        title('$K$','Interpreter', 'latex')
        legend('Correlated', 'IID','Misspecified');
    end
    if k==3
        ylim([0 12])
        xtickformat('%.1f');
        ytickformat('%.0f');
        xticks([1.0 2.0 3.0 4.0]);
        yticks([0 3 6 9 12]);
        title('$C_0$','Interpreter', 'latex')
    end
    box off;
end
print('Fig2_2','-dpng','-r300')