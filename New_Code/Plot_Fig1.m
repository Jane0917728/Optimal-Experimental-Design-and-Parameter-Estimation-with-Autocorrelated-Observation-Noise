% Drawing the MLE solution for the correlated noise, IID noise and
% misspecified noise (assummed IID but it's correlated noise)
clear;
clc;
r=0.2; K=50; C0=4.5;
phi=0.1; 
SigmaIID=sqrt(9);
SigmaOU=sqrt(SigmaIID^2*(2*phi));
% SigmaC=1;
% Sigma=SigmaC;
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
 
lb=[0.05,35,0.1];
ub=[0.5,65,15]; 
% 初始化 lower 和 upper 数组
lower = 100*ones(types,length(t1)); % 预测区间的下限
upper = zeros(size(lower)); % 预测区间的上限
M=1000;
 current_lower=zeros(M,length(t1));
 current_upper=current_lower;
M = 1000; % Number of parameter samples to be drawn
sampled = zeros(M, num_param); % Initialize array to store samples of [r, k, C0] (3 columns)

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

     
        lls = zeros(M, 1); % Initialize array to store log-likelihood values for each sample
        kount = 0; % Initialize sampling counter

        % Parameter range definition (lower and upper bounds)
        lb1 = [0.1, 45, 2];   % Lower bounds for [r, k, N0]
        ub1= [0.3, 55, 8];  % Upper bounds for [r, k, N0]
 
        % Sampling process
        while kount < M
            % Randomly draw values for r, k, and N0 from uniform distributions within bounds
            sample = lb1 + rand(1, 3) .* (ub1 - lb1); % Draw [r, k, N0] within bounds

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
      


        % 对每个样本进行处理        
            %Ts = @(t) rsampled(i) + (N0 - rsampled(i)) * exp(-ksampled(i) * t);     
            for i=1:M 
                        % 计算当前时间点的预测值加上噪声模型的上下界
                current_lower(i,:) = Nfunction(sampled(i,1),sampled(i,2),t1,sampled(i,3)) + lower_bound;
                current_upper(i,:) =Nfunction(sampled(i,1),sampled(i,2),t1,sampled(i,3))  + upper_bound;
            end
                % Update lower and upper  
                lower(k,:)=min(current_lower,[],1);
                upper(k,:)=max(current_upper,[],1);
         %function [solution,obj_value]=Opt_prediction_min(x0,lb,ub)
           fig=figure('Position',[20 20 1250 320],'color','w');        
        subplot('Position',[0.06,0.14,0.23,0.75]);
        hold on;
        fill([t1, fliplr(t1)], [lower(k,:), fliplr(upper(k,:))], [0.8, 1, 0.8], 'EdgeColor', 'none');

        h1=plot(t1,posterior(:,k), '-',t1, N_real1,'--',t,  measurements(:,k),'o', 'MarkerFaceColor', 'c', 'MarkerSize',6.5,'LineWidth', 2);
%      ylim=([min(lower(k,:))-1,max( upper(k,:))+1]);
       xlabel('time');           
        ylabel('Population');
        axis('square');           
        legend(h1,'MLE solution','real','measurement data','Location', 'southeast');
%         ylim([min(min(measurements(:,k)),min(prediction_min(:,k)))-1,...
%         max(max(measurements(:,k)),max(prediction_max(:,k)))+1]);
    
%          ylim=(h1,[min(lower(k,:))-1,max( lower(k,:))+1]);
            
        title(sprintf('%s', type{k})); 
         
        hold off;
                 
             

                
           
        
      
         %--------------------------------------------------------------------------------------------------------
      %% profile likelihood
      % An array of parameters of interest (assign the value)
      param_vals=zeros(num_param,numpts);
      % Likelihood array
       max_ls=zeros(num_param,numpts); 
       % Optimized parameter array
       minimizers=zeros(num_param,numpts); 
%         lb=[0.07,40,1];
%         ub=[0.45,60,12]; 
       optimal_param_vals=mlevalue';
       param_vals=[param_vals,optimal_param_vals];
 

        %%%for r
         initial_likelihood=optimal_param_vals;
        for param=1:num_param    
        % The parameters are equally divided and the index of the globally optimal parameter is obtained
            param_vals(param,1:numpts)=linspace(lb(param),ub(param),numpts); % assign values to a parameter.
            param_vals(param,:)=sort(param_vals(param,:));
            [~,mle_idx]=min(abs(param_vals(param,:)-optimal_param_vals(param)));% find the index of the optimal parameter
           fixed=fixed1;
           fixed(param)=1;    
           minimizers(fixed==0,mle_idx)=optimal_param_vals(fixed==0);
        %    minimizers(:,mle_idx)=[optimal_param_vals(2),optimal_param_vals(1)];
           max_ls(param,mle_idx)=-Max_l;
        %    MSE_v(param,mle_idx)=mse;

            for i=mle_idx+1:numpts+1
                  % if optimize K, given value of r, let initial K=optimal value of K in previous step 
               initial_likelihood(fixed==0)=minimizers(fixed==0,i-1); 
                      % Let r=given value defined in param_value.
                initial_likelihood(param)=param_vals(param,i);
                [minimizer,max_likelihood,~,~] =  fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements(:,k)',t);

                %     [minimizer,MSE_v(param,i),max_ls(param,i),~,~,~,~] =fmincon_likelihood2(initial_likelihood,fixed_params,lb,ub,N0,t,NData,Sigma);
                minimizers(fixed==0,i)=minimizer;
                max_ls(param,i)=-max_likelihood;
        %         mse=(Nfunction(minimizer(1),minimizer(2),t,N0)-Nfunction(r,K,t,N0)).^2/num_steps;
            end

            for i=mle_idx-1:-1:1
        %         fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
                initial_likelihood(fixed==0)=minimizers(fixed==0,i+1);
                initial_likelihood(param)=param_vals(param,i);
                [minimizer,max_likelihood,~,~] = fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements(:,k)',t);
                 minimizers(fixed==0,i)=minimizer;
                 max_ls(param,i)=-max_likelihood;
            end
        end %end of num_param


%%


        % Stores the zero of the Profile Likelihood curve
       free_param_count=1;
      % Stores the zero of the Profile Likelihood curve
      zs = cell(num_param,1);
      % Store The width of the store confidence interval
      conf_interval=nan(num_param,1);
        %  [left bottom width height] 
        %  fig=figure('Position',[20 20 1000 320],'color','w');
        
%         fig=figure('Position',[20 20 1250 320],'color','w');
%         
        for param=1:num_param    

         free_param_count = free_param_count+1;
     subplot('Position',[0.06+param*0.24,0.14,0.22,0.75]);
        %     subplot(1,num_param+1,free_param_count);
        %      subplot('Position',[0.05+0.24*param,0.14,0.23,0.8]);%   [left bottom width height]

            hold on;
            xx=param_vals(param,:);
            yy=max_ls(param,:)-max(max_ls(param,:));
            plot(xx,yy,'-','color',[0,0,0],'LineWidth',2);
            xline(mlevalue(1,param),'color',[1.0 0 0],'LineWidth',1)
            plot([min(param_vals(param,:)),max(param_vals(param,:))],[-1.92,-1.92]);    
            % plot([min(param_vals(param,:)),max(param_vals(param,:))]);
            xlabel(prm_name{param});
            %     ylabel('log(L)');
            if param==1
               ylabel('Normalized Likelihood');
            end
           axis('square');
           xlim([min(param_vals(param,:)),max(param_vals(param,:))]);
           ylim([-2.5,0]);
        %       ylim([-10,0]);
            hold off;    
            zs{param}=interp_zero(xx,yy+1.92);
            if size(zs{param},2) == 2
                 if param==1
                        CI_r(1,k)=zs{param}(1);
                        CI_r(3,k)=zs{param}(2);
                 elseif param==2
                        CI_K(1,k)=zs{param}(1);
                        CI_K(3,k)=zs{param}(2);
                 else
                        CI_N0(1,k)=zs{param}(1);
                        CI_N0(3,k)=zs{param}(2);
                 end             
                conf_interval(param)=zs{param}(2)-zs{param}(1);
                fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',prm_name{param});
                fprintf('width=%.4f: [%.4f,%.4f]\n',conf_interval(param),zs{param}(1),zs{param}(2));
                xline(zs{param}(1), 'color',[0 1 0]);  
                xline(zs{param}(2), 'color',[0 1 1]); 
            else
                fprintf('Do not have 2 intercepts for param %s, they are:\n',prm_name{param});
                disp(zs{param});
                if size(zs{param},2) == 1
                    if(zs{param} > mlevalue(1,param))
                        xline(zs{param}(1), 'color',[0 1 1]);  
                    else
                        xline(zs{param}(1), 'color',[0 1 0]);  
                    end
                end
                %%-------------------------------------------------------------------------------------------------------------begin
                 if size(zs{param},2) == 1
                        if(zs{param} > mlevalue(1,param))
                            if param==1
                                CI_r(1,k)=-0.2;
                                CI_r(3,k)=zs{param}(1);
                            elseif param==2
                                CI_K(1,k)=10;
                                CI_K(3,k)=zs{param}(1);
                            else
                                CI_N0(1,k)=-1;
                                CI_N0(3,k)=zs{param}(1);
                            end    
                        else
                            if param==1
                                CI_r(1,k)=zs{param}(1);
                                CI_r(3,k)=1;
                            elseif param==2
                                CI_K(1,k)=zs{param}(1);
                                CI_K(3,k)=110;
                            else
                                CI_N0(1,k)=zs{param}(1);
                                CI_N0(3,k)=30;
                            end    
                        end

                    else 
                        if param==1
                                CI_r(1,k)=-0.2;
                                CI_r(3,k)=1;
                        elseif param==2
                                CI_K(1,k)=10;
                                CI_K(3,k)=110;
                        else
                                CI_N0(1,k)=-1;
                                CI_N0(3,k)=30;
                         end            
                 end
                 %------------------------------------------------------end
 
                 
                 
             end
        end
         
%--------------------------------------------------------------------------------------------------------------------------
% draw the prediction intervals

  
  

 end
  
 

 


        
      
 
   
