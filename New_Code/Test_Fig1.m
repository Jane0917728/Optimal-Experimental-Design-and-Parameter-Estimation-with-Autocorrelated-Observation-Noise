% Drawing the MLE solution for the correlated noise, IID noise and
% misspecified noise (assummed IID but it's correlated noise)
clear;
clc;
r=0.2; K=50; C0=4.5;
phi=0.02; 
SigmaIID=sqrt(9);
SigmaOU=sqrt(SigmaIID^2*(2*phi));
t_initial=0;
TF=80+t_initial;

%Estimate the parameters 'prm_interest' from 'prm_name_all'.
prm_all=[r K C0]'; prm_name_all = {'r', 'K', 'C_0'};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

numpts=70; 

n_s=11;%number of measurement points
t=linspace(t_initial,TF,n_s);
draw_n=200;
t1=linspace(t_initial,TF,draw_n);
N_real=Nfunction(r,K,t,C0);
N_real1=Nfunction(r,K,t1,C0);
est_matrix=zeros(3,3);
 

SIG_OU=zeros(n_s,n_s); % The covariance matrix 
        for s = 1:n_s
            for tindex = 1:n_s
%                 SIG(s, tt) = Sigma^2*(exp(-phi*abs(t_opt(tt)-t_opt(s)))-exp(-phi*(t_opt(tt)+t_opt(s))))/(2*phi);   
               SIG_OU(s, tindex) = SigmaOU^2*exp(-phi*abs(t(tindex)-t(s)))/(2*phi);   
            end   
        end 
 SIG_IID=eye(size(SIG_OU))*SigmaIID^2;  
 types=3;
 type={'OU noise','IID noise','Misspecified noise'};
 Z=randn(1,n_s);
%  Z=[0.825218894	1.378971978	-1.058180258	-0.468615581	-0.272469409	1.098424618	-0.277871933	0.701541458	-2.0518163	-0.353849998	-0.823586525
% ];
 OU_noise=Generate_OU_seeds(phi,SigmaOU,t,Z);
 IID_noise=Z(1:end)*SigmaIID;
 initial=[0.1,60,10]; 
 posterior=zeros(draw_n,types);
%  define_of_CI
CI_r=zeros(3,types);
CI_K=CI_r;
CI_N0=CI_r;


measurements=zeros(n_s,types);
 lower_bound = norminv(0.05, 0, SigmaIID); %use for predication interval
upper_bound = norminv(0.95, 0, SigmaIID); %use for predication interval
 
lb=[0.05,35,0.5];
ub=[0.45,65,12]; 
% 初始化 lower 和 upper 数组
lower = 100*ones(length(t1),types); % 预测区间的下限
upper = zeros(size(lower)); % 预测区间的上限


 for k=1:types    
      fixed1=[0,0,0]; fixed=fixed1;
      if k==1
         %OU noise
          SIG=SIG_OU;   
          noise=OU_noise;       
      else
          if k==2
              SIG=SIG_IID;
              noise=IID_noise;
          else
              SIG=SIG_IID;
              noise=OU_noise;
          end
      end  
       measurements(:,k)=N_real'+noise';
      
       [mlevalue,max_likelihood] =fmincon_likenew(fixed,initial,lb,ub,SIG,measurements(:,k)',t);
    %    prm_all=mlevalue;
       est_matrix(k,:)=mlevalue;
       CI_r(2,k)=mlevalue(1);
       CI_K(2,k)=mlevalue(2);
       CI_N0(2,k)=mlevalue(3);
 
      posterior(:,k)=Nfunction(mlevalue(1),mlevalue(2),t1',mlevalue(3));
      Max_l=max_likelihood;
      
      
      % Sampling the parameters
        % Parameter definition
        df = 3; % Degrees of freedom set to 3 (assuming degrees of freedom are directly related to the number of parameters)
        llstar = -chi2inv(0.95, df) / 2; % Log-likelihood threshold

        M = 1000; % Number of samples to be drawn
        sampled = zeros(M, 3); % Initialize array to store samples of [r, k, N0] (3 columns)
        lls = zeros(M, 1); % Initialize array to store log-likelihood values for each sample
        kount = 0; % Initialize sampling counter

        % Parameter range definition (lower and upper bounds)
%         lb = [0.1, 40, 1];   % Lower bounds for [r, k, N0]
%         ub = [0.3, 60, 10];  % Upper bounds for [r, k, N0]

        % Sampling process
        while kount < M
            % Randomly draw values for r, k, and N0 from uniform distributions within bounds
            sample = lb + rand(1, 3) .* (ub - lb); % Draw [r, k, N0] within bounds

            % Calculate the log-likelihood for the current sample and subtract fmle (assumed known or calculated through other means)
%             current_ll = loglhood(data, sample(3), sample(1:2), sigma, t) - fmle;
              current_ll=- est_likenew(sample,SIG,measurements(:,k)',t) - (-Max_l);
     
            % If the current sample meets the condition, retain the sample
            if current_ll >= llstar
                kount = kount + 1; % Update the counter

                % Store the sample values in the sampled array
                sampled(kount, :) = sample; % Store [r, k, N0]
                lls(kount) = current_ll;    % Store log-likelihood value
            end
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      


 
        
            for j = 1:length(t1)
                for i = 1:M
                          % The prediction value + noise bounds计算当前时间点的预测值加上噪声模型的上下界
                current_lower = Nfunction(r,K,t1(j),C0) + lower_bound;
                current_upper =Nfunction(r,K,t1(j),C0) + upper_bound;

                % Update lower and upper  
                if current_lower < lower(j,k)
                    lower(j,k) = current_lower;
                end

                if current_upper > upper(j,k)
                    upper(j,k) = current_upper;
                end
            end
        end
      
         %--------------------------------------------------------------------------------------------------------
      
         
%--------------------------------------------------------------------------------------------------------------------------
% draw the prediction intervals

  
  figure;
%function [solution,obj_value]=Opt_prediction_min(x0,lb,ub)
        subplot('Position',[0.06,0.14,0.23,0.75]);
        hold on;
        fill([t1, fliplr(t1)], [lower(:,k)', fliplr(upper(:,k)')], [0.8, 1, 0.8], 'EdgeColor', 'none');

        h1=plot(t1,posterior(:,k), '-',t1, N_real1,'--',t,  measurements(:,k),'o', 'MarkerFaceColor', 'c', 'MarkerSize',6.5,'LineWidth', 2);
    
   xlabel('time');           
    ylabel('Population');
    axis('square');           
    legend(h1,'MLE solution','real','measurement data','Location', 'southeast');

    title(sprintf('%s', type{k})); 
    hold off;
             
 end